import json

claims = [json.loads(l) for l in open("results/claims.jsonl")]
active = [c for c in claims if not c.get("retracted")]
retracted = [c for c in claims if c.get("retracted")]
iter9 = [c for c in active if c.get("iteration") == 9]

print(f"Total: {len(claims)}, Active: {len(active)}, Retracted: {len(retracted)}, Iter9: {len(iter9)}")
for c in iter9:
    ct = c["claim_type"]
    expr = c["expression"][:90]
    print(f"  [{ct}] {expr}")

gov = [json.loads(l) for l in open("results/governance.jsonl")]
print(f"\nGovernance records: {len(gov)}")
for g in gov:
    evt = g.get("event", g.get("mode", "?"))
    mode = g.get("mode", "?")
    it = g.get("iteration", g.get("iteration_retracted", "?"))
    print(f"  iter={it} event={evt} mode={mode}")
