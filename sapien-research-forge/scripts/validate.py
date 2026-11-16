import json, pathlib, unittest, math
r=pathlib.Path(__file__).resolve().parents[1]
d=json.loads((r/"registry/results.json").read_text(encoding="utf-8"))
src=json.loads((r/"evidence/source_manifest.json").read_text(encoding="utf-8"))
assert len(set(x["result_id"] for x in d["results"]))==len(d["results"])
for item in src["files"]:
    assert (r/item["path"].removeprefix("sapien-research-forge/")).is_file(),item["path"]
for phi in (0.1,0.9,1.7):
    eps=0.005
    actual=sum((1+eps*math.cos(phi-2*math.pi*k/3))**(-0.5) for k in range(3))
    approx=3+9*eps**2/16-15*eps**3*math.cos(3*phi)/64+315*eps**4/1024
    assert abs(actual-approx)<2*eps**5
print("PASS",len(d["results"]),"claims,",len(src["files"]),"source transcriptions, D3 smoke tests.")
