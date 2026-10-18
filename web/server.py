#!/usr/bin/env python3
"""Live dashboard for smol runs.  python3 server.py [--runs ~/smol/runs] [--port 8095] [--host 0.0.0.0]

  /                 dashboard: every condition and run, refreshed every 30 s
  /api/runs         JSON: status and population series per run (from log.txt)
  /viewer?run=C/S   the grid viewer for one run, built on demand and cached

Read-only. Standard library only.
"""
import argparse, html, json, os, re, subprocess, sys, threading, time
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--runs", default=os.path.expanduser("~/smol/runs"))
ap.add_argument("--port", type=int, default=8095)
ap.add_argument("--host", default="0.0.0.0")
ap.add_argument("--viz", default=os.path.join(HERE, "..", "viz", "build.py"))
ap.add_argument("--cache", default=os.path.expanduser("~/smol/viewcache"))
A = ap.parse_args()
os.makedirs(A.cache, exist_ok=True)

STEP = re.compile(r"^Step (\d+) done, number of molecules is (\d+), nbound = (\d+), (\d+) steps/s(?: \| (.*))?")
EXT = re.compile(r"^EXTINCT at step (\d+)")
SWITCHES = ("DISPLACE", "DISPLACE_RATE", "DISPLACE_K", "BINDLEN", "TWOSIDED", "COMPAT_PARTNER",
            "COMPAT_DECAY", "DECAY_INDEPENDENT", "ENCOUNTER_TURN", "DECAY", "MUTATE", "NSTEPS")
_cache = {}
_lock = threading.Lock()


def parse_run(path):
    log = os.path.join(path, "log.txt")
    if not os.path.exists(log):
        return {"status": "pending", "series": [], "step": 0, "n": 0, "bound": 0}
    st = os.stat(log)
    key = (log, st.st_mtime, st.st_size)
    with _lock:
        if _cache.get(log, (None,))[0] == key:
            return _cache[log][1]
    series, last, ext, counters = [], None, None, {}
    with open(log, errors="replace") as f:
        for line in f:
            m = STEP.match(line)
            if m:
                s, n, b, rate = map(int, m.groups()[:4])
                series.append((s, n, b))
                last = (s, n, b, rate)
                if m.group(5):
                    p = m.group(5).split()
                    counters = {p[i]: int(p[i + 1]) for i in range(0, len(p) - 1, 2)}
            else:
                e = EXT.match(line)
                if e:
                    ext = int(e.group(1))
    conf = {}
    rc = os.path.join(path, "run.conf")
    if os.path.exists(rc):
        for line in open(rc):
            p = line.split()
            if len(p) >= 2 and p[0] in SWITCHES:
                conf[p[0]] = p[1]
    nsteps = int(float(conf.get("NSTEPS", 1000000)))
    age = time.time() - st.st_mtime
    if ext is not None:
        status = "extinct"
    elif last and last[0] >= nsteps:
        status = "done"
    elif age < 120:
        status = "running"
    else:
        status = "stopped"
    # downsample to <= 250 points
    k = max(1, len(series) // 250)
    ds = series[::k]
    if series and ds[-1] != series[-1]:
        ds.append(series[-1])
    r = {"status": status, "step": ext if ext is not None else (last[0] if last else 0),
         "n": last[1] if last else 0, "bound": last[2] if last else 0, "rate": last[3] if last else 0,
         "series": ds, "counters": counters, "conf": conf, "nsteps": nsteps, "age": int(age)}
    with _lock:
        _cache[log] = (key, r)
    return r


def descriptions():
    p = os.path.join(A.runs, "conditions.json")
    try:
        return json.load(open(p))
    except Exception:
        return {}


def all_runs():
    out = []
    desc = descriptions()
    for cond in sorted(os.listdir(A.runs)):
        cdir = os.path.join(A.runs, cond)
        if not os.path.isdir(cdir):
            continue
        for seed in sorted(os.listdir(cdir), key=lambda s: int(re.sub(r"\D", "", s) or 0)):
            p = os.path.join(cdir, seed)
            if os.path.isdir(p):
                r = parse_run(p)
                r.update(cond=cond, seed=seed, desc=desc.get(cond, ""))
                out.append(r)
    return out


def build_viewer(run):
    cond, seed = run.split("/", 1)
    rdir = os.path.join(A.runs, cond, seed)
    if not (os.path.isdir(rdir) and os.path.exists(os.path.join(rdir, "grid.snap"))):
        return None, "no such run, or it has no grid.snap yet"
    out = os.path.join(A.cache, f"{cond}__{seed}.html")
    snap = os.path.join(rdir, "grid.snap")
    fresh = os.path.exists(out) and os.path.getmtime(out) >= os.path.getmtime(snap) - 600
    if not fresh:
        r = parse_run(rdir)
        every = max(500, (r["step"] // 400) // 500 * 500)
        p = subprocess.run([sys.executable, A.viz, rdir, out + ".tmp", str(every)], capture_output=True, text=True)
        if p.returncode != 0:
            return None, p.stderr[-2000:]
        os.replace(out + ".tmp", out)
    return out, None


class H(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def send(self, code, body, ctype="text/html; charset=utf-8"):
        b = body.encode() if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(b)

    def do_GET(self):
        u = urlparse(self.path)
        if u.path == "/":
            self.send(200, open(os.path.join(HERE, "dashboard.html"), "rb").read())
        elif u.path == "/api/runs":
            self.send(200, json.dumps(all_runs()), "application/json")
        elif u.path == "/viewer":
            run = parse_qs(u.query).get("run", [""])[0]
            if ".." in run or run.count("/") != 1:
                return self.send(400, "bad run")
            path, err = build_viewer(run)
            if err:
                return self.send(404, f"<pre>{html.escape(err)}</pre>")
            self.send(200, open(path, "rb").read())
        else:
            self.send(404, "not found")


print(f"serving {A.runs} on http://{A.host}:{A.port}/", flush=True)
ThreadingHTTPServer((A.host, A.port), H).serve_forever()
