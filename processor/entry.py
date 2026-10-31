"""
One image, two runtimes. Fixes the thread count from THREADS (safe_inputs.cpu_threads: else the
cgroup CPU quota, else the CPUs this process may use) before numpy is imported; then, as a Lambda
function (AWS_LAMBDA_RUNTIME_API is set), hands over to the Lambda runtime interface client with
processor.handler, or, as a container, runs main() once.
"""
import os
import sys

from processor.safe_inputs import cpu_threads

raw = os.environ.get("THREADS", "").strip()
n = cpu_threads(int(raw) if raw.isdigit() else 0)
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[var] = str(n)
os.environ["SC_MARKERS_THREADS"] = str(n)

if os.environ.get("AWS_LAMBDA_RUNTIME_API"):
    os.execvp(sys.executable, [sys.executable, "-m", "awslambdaric", "processor.handler.handler"])

from processor.main import main  # noqa: E402

sys.exit(main())
