import pyarrow.feather as pf
import sys, os
for f in sys.argv[1:]:
    print("reading", f, flush=True)
    t = pf.read_table(f)
    out = f.replace(".feather", ".v2.feather")
    print("writing", out, t.num_rows, "rows", flush=True)
    pf.write_feather(t, out, version=2)
    print("done", os.path.getsize(out), flush=True)
