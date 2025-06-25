#!/bin/bash
#OAR -l nodes=4,walltime=2:00:00
#OAR -p cluster='grisou'
set -euxo pipefail

module load python/3.10
export PATH="$HOME/.local/bin:$PATH"

PORT=6379
HEAD_IP=$(hostname -i)
LOG_DIR="$HOME/ray_logs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"

function start_head() {
  echo "[HEAD] Démarrage sur $HEAD_IP:$PORT"
  ray start --head --port=$PORT --num-cpus=$(nproc) \
    > "$LOG_DIR/head.log" 2>&1
}

function start_workers() {
  for node in $(sort -u "$OAR_NODEFILE" | grep -v "$(hostname -s)"); do
    echo "[WORKER] Démarrage sur $node"
    oarsh "$node" "module load python/3.10 && \
      export PATH=$HOME/.local/bin:\$PATH && \
      ray start --address=$HEAD_IP:$PORT --num-cpus=\$(nproc) \
      > $LOG_DIR/worker_\$(hostname -s).log 2>&1"
  done
}

start_head
sleep 10
start_workers
sleep 30

echo "[TEST] Vérification du cluster"
python3.10 - <<EOF
import ray
ray.init(address="auto")
print("Nodes:", ray.nodes())
print("Ressources:", ray.cluster_resources())
ray.shutdown()
EOF

exec python3.10 kefsi_mekhazni_workspace/src/STRING_levure.py
