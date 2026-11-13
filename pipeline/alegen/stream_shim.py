"""Local OpenAI-compatible shim that turns non-streaming chat completions into streamed upstream calls.

Why it exists: on this build host every outbound connection goes through a relay that closes a connection
after about 30 seconds without response bytes. Harbor's terminus-2 agent makes non-streaming calls, and a
reasoning model often takes longer than that before the first byte, so the call dies with HTTP 502.
The shim listens on localhost (not proxied), forwards each request upstream with `stream: true`, collects
the chunks and answers with an ordinary completion object. On a host with normal networking it is not needed.

  python stream_shim.py <port> <upstream_base> [extra_json]
      e.g.  python stream_shim.py 8765 https://api.moonshot.ai/v1
            python stream_shim.py 8770 https://openrouter.ai/api/v1 '{"provider":{"only":["moonshotai"],"allow_fallbacks":false}}'

`extra_json` is merged into every request body. It exists for aggregators: OpenRouter serves one model name from
many providers, several of them with quantised weights, and picks one per request unless the request pins it.

The upstream credential is not handled here: the Authorization header of the caller is passed on unchanged
(on this host the network proxy replaces it with the saved credential).
One line per request is appended to <port>.shim.log next to the current directory (status, seconds, tokens).
"""
from __future__ import annotations

import json
import sys
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import httpx

UPSTREAM = ""
EXTRA: dict = {}
LOG = None


def log(msg: str) -> None:
    line = time.strftime("%H:%M:%S ") + msg
    if LOG:
        LOG.write(line + "\n")
        LOG.flush()


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, *a):  # silence the default access log
        pass

    def _send(self, status: int, body: bytes, ctype: str = "application/json") -> None:
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _headers(self) -> dict:
        h = {"Content-Type": "application/json"}
        if self.headers.get("Authorization"):
            h["Authorization"] = self.headers["Authorization"]
        return h

    def do_GET(self):
        try:
            r = httpx.get(UPSTREAM + self.path.removeprefix("/v1"), headers=self._headers(), timeout=60)
            self._send(r.status_code, r.content)
        except Exception as e:
            self._send(502, json.dumps({"error": {"message": f"shim: {e}"}}).encode())

    def do_POST(self):
        n = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(n)
        t0 = time.time()
        try:
            req = json.loads(raw)
        except Exception:
            self._send(400, b'{"error":{"message":"shim: body is not JSON"}}')
            return
        if not self.path.endswith("/chat/completions") or req.get("stream"):
            self._send(400, b'{"error":{"message":"shim: only non-streaming chat completions are supported"}}')
            return
        # Some providers reject a history that contains an assistant turn without text (the agent stores such a
        # turn when the model returned nothing visible). Give those turns a visible marker instead of failing the trial.
        for m in req.get("messages") or []:
            if m.get("role") == "assistant" and not m.get("tool_calls") and not (m.get("content") or "").strip():
                m["content"] = "(empty response)"
        req.update(EXTRA)
        req["stream"] = True
        req["stream_options"] = {"include_usage": True}
        # Transient upstream failures (overload, gateway errors, a stream cut before the end) are retried here with
        # a growing pause, so that one busy minute at the provider does not end a two-hour trial.
        for attempt, pause in enumerate((0, 5, 15, 30, 60, 90, 120)):
            if pause:
                time.sleep(pause)
            content, reasoning, finish, usage, meta = [], [], None, None, {}
            status, body, err = 200, b"", None
            try:
                with httpx.stream("POST", UPSTREAM + "/chat/completions", headers=self._headers(), json=req,
                                  timeout=httpx.Timeout(connect=30, read=120, write=60, pool=30)) as r:
                    status = r.status_code
                    if status != 200:
                        body = r.read()
                    else:
                        for line in r.iter_lines():
                            if not line or not line.startswith("data:"):
                                continue
                            data = line[5:].strip()
                            if data == "[DONE]":
                                break
                            ch = json.loads(data)
                            for k in ("id", "model", "created", "system_fingerprint", "provider"):
                                if ch.get(k) is not None:
                                    meta[k] = ch[k]
                            if ch.get("usage"):
                                usage = ch["usage"]
                            for c in ch.get("choices") or []:
                                d = c.get("delta") or {}
                                if d.get("content"):
                                    content.append(d["content"])
                                if d.get("reasoning_content"):
                                    reasoning.append(d["reasoning_content"])
                                if c.get("finish_reason"):
                                    finish = c["finish_reason"]
                                if c.get("usage"):
                                    usage = c["usage"]
            except Exception as e:
                err = f"{type(e).__name__}: {e}"
            if status == 200 and err is None and finish is not None:
                break
            what = err or (f"HTTP {status}: {body[:160]!r}" if status != 200 else "stream ended without finish_reason")
            # an exhausted account also answers 429 (or 402); waiting does not help, and retrying hid it for over an
            # hour once. A rate limit is different: Google words it "exceeded your current quota ... check your plan
            # and billing details" and lifts it within a minute, so "billing" alone must not count as exhausted.
            fatal = status == 402 or (status == 429 and any(k in body.lower() for k in (b"insufficient balance", b"suspended")))
            retryable = not fatal and (err is not None or status in (408, 409, 429, 500, 502, 503, 504, 529)
                                       or (status == 200 and finish is None))
            log(f"attempt {attempt + 1} failed after {time.time() - t0:.1f}s: {what}" + ("" if retryable else " (not retried)"))
            if not retryable:
                self._send(status, body)
                return
        else:
            if status != 200 and body:
                self._send(status, body)
            else:
                self._send(502, json.dumps({"error": {"message": f"shim: upstream failed after retries: {what}"}}).encode())
            return
        msg = {"role": "assistant", "content": "".join(content)}
        if reasoning:
            msg["reasoning_content"] = "".join(reasoning)
        out = {"id": meta.get("id", "shim"), "object": "chat.completion", "created": meta.get("created", int(t0)),
               "model": meta.get("model", req.get("model")),
               "choices": [{"index": 0, "message": msg, "finish_reason": finish}],
               "usage": usage or {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0}}
        cost = f" cost={usage['cost']}" if usage and usage.get("cost") is not None else ""   # billed amount, if the upstream reports it
        prov = f" provider={meta['provider']!r}" if meta.get("provider") else ""
        log(f"ok {time.time() - t0:.1f}s finish={finish}{prov}{cost} usage={json.dumps(usage)[:160] if usage else None}")
        self._send(200, json.dumps(out).encode())


def main() -> int:
    global UPSTREAM, EXTRA, LOG
    port, UPSTREAM = int(sys.argv[1]), sys.argv[2].rstrip("/")
    EXTRA = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    LOG = open(f"{port}.shim.log", "a")
    log(f"listening on 127.0.0.1:{port} -> {UPSTREAM}" + (f" extra={json.dumps(EXTRA)}" if EXTRA else ""))
    ThreadingHTTPServer(("127.0.0.1", port), Handler).serve_forever()
    return 0


if __name__ == "__main__":
    sys.exit(main())
