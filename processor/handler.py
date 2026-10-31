"""
The Lambda runtime. A lambda step's payload carries inputDir, outputDir and
the step's parameters as fields (where a container gets them as environment
variables); the handler sets those variables and runs the same main(), so
both runtimes write the same outputs. The thread count was fixed by
processor/entry.py when the function started (THREADS, else 2: Lambda exposes
no cgroup CPU quota; 3008 MB gets about two vCPUs).
"""
import os

from processor import main as app


def handler(event, context):
    os.environ["INPUT_DIR"] = event.get("inputDir") or ""
    os.environ["OUTPUT_DIR"] = event.get("outputDir") or ""
    # A warm function keeps its environment between invocations, so a
    # parameter this run doesn't set is cleared, not left from the last one.
    for name in app.PARAMETERS:
        value = event.get(name)
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = str(value).lower() if isinstance(value, bool) else str(value)
    code = app.main()
    if code != 0:
        # An exception, so the step fails; a returned value would succeed.
        raise RuntimeError(f"sc-markers exited with {code}")
    return {"status": "ok", "outputDir": os.environ["OUTPUT_DIR"]}
