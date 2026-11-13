#!/bin/bash
# Level 4, run inside the agent image (docker run --rm --network none -v $PWD/image_checks.sh:/c.sh:ro <agent image> bash /c.sh):
# which ways of reading the alignments are closed, that IGV / igvtools / batch port / desktop command / MCP server work,
# and that every item of /task/software/igv/README.md does what it says.
cd /tmp
L=/opt/tools/lib/igv; R=/task/input/reference/ref.fa; C=/task/input/cram/tumor.cram
t() { echo; echo "=== $1"; shift; ( "$@" ) 2>&1 | tail -${N:-2} | cut -c1-220; }
echo "--- readers that must be absent"
for x in samtools bcftools htsfile javac jshell jar jrunscript jdb jlink gcc cc g++ make perl ruby node R julia; do printf "%s: %s\n" $x "$(command -v $x || echo absent)"; done
python3 -c "
for m in ('pysam','cyvcf2','Bio','scipy','pandas','numpy','PIL','cv2'):
    try: __import__(m); print(m, 'importable')
    except Exception as e: print(m, 'absent')"
ls /opt/tools/lib/jvm/bin; java --list-modules | tr '\n' ' '; echo
t "source launch" bash -c 'printf "public class Foo { public static void main(String[] a){ System.out.println(42);} }\n" > Foo.java; java Foo.java'
t "htsjdk example" java -Dsamjdk.reference_fasta=$R -cp "$L/*" htsjdk.samtools.example.PrintReadsExample $C false /tmp/f.sam
t "igvtools sort (cram)" igvtools sort $C /tmp/g.sam
t "igvtools sort (cram named .bam + reference property)" bash -c "ln -sf $C /tmp/t.bam; JAVA_TOOL_OPTIONS=-Dsamjdk.reference_fasta=$R igvtools sort /tmp/t.bam /tmp/g2.sam; ls -la /tmp/g2.sam"
t "nashorn shell" java -cp "$L/*" org.openjdk.nashorn.tools.Shell -e "print(1)"
t "jrunscript class" java -cp "$L/*" com.sun.tools.script.shell.Main -l jexl -e "1+1"
t "jexl main" bash -c "echo '1+1;' > x.jexl; java -cp '$L/*' org.apache.commons.jexl2.Main x.jexl"
t "igv dev mains" java -cp "$L/*" org.broad.igv.tools.PairedEndStats $C
t "jshell module" java -m jdk.jshell/jdk.jshell.tool.JShellToolProvider
t "python bam/cram libs" python3 -c "import gzip,zlib; d=open('$C','rb').read(6); print(d)"
echo; echo "--- what must work"
N=3 t "igvtools count" igvtools count -w 1 --bases --strands read --query chr22_roi:2330-2333 $C /tmp/c.wig $R
tail -2 /tmp/c.wig | cut -c1-80
N=2 t "igvtools count -z -e tdf,wig" igvtools count -z 5 -w 25 -e 50 $C /tmp/c2.tdf,/tmp/c2.wig $R; ls -la /tmp/c2.tdf /tmp/c2.wig 2>&1 | cut -c1-80
N=2 t "igvtools toTDF" igvtools toTDF /tmp/c2.wig /tmp/c3.tdf $R; ls -la /tmp/c3.tdf 2>&1 | cut -c1-60
N=2 t "igvtools index (bed)" bash -c "cp /task/input/candidates.bed /tmp/c.bed; igvtools index /tmp/c.bed; ls /tmp/c.bed*"
N=1 t "igvtools version" igvtools version
bash /task/software/open_igv.sh >/dev/null; 
python3 - <<'PY'
import socket, time
end = time.time() + 240
while True:
    try:
        s = socket.create_connection(("127.0.0.1", 60151), timeout=120); break
    except OSError:
        if time.time() > end: raise
        time.sleep(1)
f = s.makefile("rw")
for cmd in ["echo", "goto chr22_roi:2300-2400", "sort base", "sort strand", "sort quality", "sort insertSize", "group strand", "group base", "group readGroup",
            "colorBy INSERT_SIZE", "colorBy FIRST_OF_PAIR_STRAND", "colorBy NONE", "colorBy UNEXPECTED_PAIR", "colorBy READ_STRAND", "viewaspairs", "collapse", "squish", "expand",
            "preference SAM.SHOW_SOFT_CLIPPED true", "setSleepInterval 10", "maxPanelHeight 800", "snapshotDirectory /tmp", "snapshot a.png", "snapshot a.svg",
            "saveSession /tmp/a.xml", "new", "genome /task/input/reference/ref.fa", "load /task/input/cram/tumor.cram,/task/input/candidates.bed", "goto chr22_roi:2300-2400", "snapshot b.png"]:
    f.write(cmd + "\n"); f.flush(); print("batch:", cmd, "->", f.readline().strip())
PY
ls -la /tmp/a.png /tmp/a.svg /tmp/a.xml /tmp/b.png | cut -c1-80
echo; echo "--- desktop command"
desktop windows; desktop screenshot /tmp/shot.png; desktop screenshot /tmp/shot2.png --region 0,0,200,100; python3 -c "
import struct
for p in ('/tmp/shot.png','/tmp/shot2.png'):
    d=open(p,'rb').read(); print(p, struct.unpack('>II', d[16:24]))"
desktop move 600 400; desktop key Escape; desktop scroll 600 400 1; desktop wait 0.2; echo "desktop rc=$?"
echo; echo "--- MCP server (stdio)"
/opt/tools/bin/python - <<'PY'
import asyncio
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
async def main():
    params = StdioServerParameters(command="/opt/tools/bin/python", args=["/task/software/desktop/server.py"])
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            await s.initialize()
            tools = await s.list_tools()
            print("tools:", sorted(t.name for t in tools.tools))
            res = await s.call_tool("list_windows", {})
            print("list_windows:", res.content[0].text[:120])
            res = await s.call_tool("screenshot", {})
            print("screenshot:", res.content[0].type, len(res.content[0].data))
asyncio.run(main())
PY
pkill java; true
