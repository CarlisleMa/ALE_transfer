# Lock-down of the review workstation image: IGV keeps working, everything else that could read
# alignments through the Java runtime or IGV's bundled libraries is removed.
import glob, os, shutil, struct, subprocess, zipfile
T = "/opt/tools"
JVM, LIB = T + "/lib/jvm", T + "/lib/igv"
MODS = ("java.base,java.datatransfer,java.desktop,java.instrument,java.logging,java.management,java.naming,"
        "java.net.http,java.prefs,java.security.jgss,java.security.sasl,java.sql,java.xml,java.xml.crypto,"
        "jdk.charsets,jdk.crypto.ec,jdk.crypto.cryptoki,jdk.localedata,jdk.management,jdk.naming.dns,jdk.net,"
        "jdk.security.auth,jdk.unsupported,jdk.unsupported.desktop,jdk.xml.dom,jdk.zipfs")
KEEP_MAIN = {"org/broad/igv/ui/Main.class", "org/broad/igv/tools/IgvTools.class"}
DROP = ("org/broad/igv/tools/sort/BAMSorter", "org/broad/igv/tools/sort/SAMSorter",
        "htsjdk/samtools/example/", "htsjdk/variant/example/", "htsjdk/tribble/example/", "htsjdk/samtools/apps/",
        "org/apache/commons/jexl2/Main", "org/apache/commons/jexl2/scripting/Main")


def demain(data):
    """Rename `public static void main(String[])` of a class file by pointing its name at another
    constant ('Code'), so the class can no longer be used as a program entry point."""
    if data[:4] != b"\xca\xfe\xba\xbe":
        return data
    n = struct.unpack(">H", data[8:10])[0]
    off, i, utf = 10, 1, {}
    while i < n:
        t = data[off]
        if t == 1:
            ln = struct.unpack(">H", data[off + 1:off + 3])[0]
            utf[i] = data[off + 3:off + 3 + ln]
            off += 3 + ln
        elif t in (3, 4, 9, 10, 11, 12, 17, 18):
            off += 5
        elif t in (5, 6):
            off += 9
            i += 1
        elif t in (7, 8, 16, 19, 20):
            off += 3
        elif t == 15:
            off += 4
        else:
            raise ValueError("constant tag %d" % t)
        i += 1
    names = {v: k for k, v in utf.items()}
    if b"main" not in names or b"([Ljava/lang/String;)V" not in names or b"Code" not in names:
        return data
    off += 6
    nif = struct.unpack(">H", data[off:off + 2])[0]
    off += 2 + 2 * nif

    def skip_members(off, patch):
        cnt = struct.unpack(">H", data[off:off + 2])[0]
        off += 2
        hits = []
        for _ in range(cnt):
            acc, name, desc, na = struct.unpack(">HHHH", data[off:off + 8])
            if patch and utf.get(name) == b"main" and utf.get(desc) == b"([Ljava/lang/String;)V":
                hits.append(off + 2)
            off += 8
            for _ in range(na):
                ln = struct.unpack(">I", data[off + 2:off + 6])[0]
                off += 6 + ln
        return off, hits
    off, _ = skip_members(off, False)
    off, hits = skip_members(off, True)
    out = bytearray(data)
    for h in hits:
        out[h:h + 2] = struct.pack(">H", names[b"Code"])
    return bytes(out)


def main():
    # 1. minimal Java runtime: no compiler, no jshell, no scripting API, no JDK tools
    subprocess.run([JVM + "/bin/jlink", "--module-path", JVM + "/jmods", "--add-modules", MODS,
                    "--no-header-files", "--no-man-pages", "--output", "/opt/igv-jre"], check=True)
    shutil.rmtree(JVM)
    shutil.move("/opt/igv-jre", JVM)
    for p in glob.glob(T + "/bin/*"):
        if os.path.islink(p) and not os.path.exists(p):
            os.unlink(p)
    # 2. IGV's bundled libraries: drop the JavaScript engine, the alignment-file sorter of igvtools
    #    and every program entry point except IGV and igvtools
    for j in glob.glob(LIB + "/nashorn-core-*.jar"):
        os.unlink(j)
    npatched = 0
    for j in sorted(glob.glob(LIB + "/*.jar")):
        src = zipfile.ZipFile(j)
        tmp = j + ".tmp"
        changed = False
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as dst:
            for info in src.infolist():
                data = src.read(info.filename)
                if info.filename.startswith(DROP):
                    changed = True
                    continue
                if info.filename.endswith(".class") and info.filename not in KEEP_MAIN:
                    new = demain(data)
                    if new != data:
                        npatched += 1
                        changed = True
                        data = new
                dst.writestr(info, data)
        src.close()
        if changed:
            os.replace(tmp, j)
        else:
            os.unlink(tmp)
    print("entry points removed:", npatched)


main()
