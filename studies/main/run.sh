
#!/usr/bin/env bash
set -Eeuo pipefail

# ANSI colours
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

# Run from the repository root, regardless of invocation directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
RUNS_DIR="$SCRIPT_DIR/runs"

cd "$REPO_ROOT"

if [[ ! -d "$RUNS_DIR" ]]; then
    printf "${RED}Error:${RESET} runs directory not found: %s\n" "$RUNS_DIR" >&2
    exit 1
fi

shopt -s nullglob
INPUTS=("$RUNS_DIR"/*/in.lammps)

if (( ${#INPUTS[@]} == 0 )); then
    printf "${RED}Error:${RESET} no in.lammps files found under %s\n" "$RUNS_DIR" >&2
    exit 1
fi

TOTAL=${#INPUTS[@]}
FAILED=0
COUNT=0
BATCH_START=$SECONDS

format_time() {
    local elapsed=$1
    printf '%02d:%02d:%02d' \
        "$((elapsed / 3600))" \
        "$(((elapsed % 3600) / 60))" \
        "$((elapsed % 60))"
}

printf "\n${BOLD}${BLUE}========================================${RESET}\n"
printf "${BOLD}${BLUE} LAMMPS batch: %d calculations${RESET}\n" "$TOTAL"
printf "${BOLD}${BLUE}========================================${RESET}\n"

for INPUT in "${INPUTS[@]}"; do
    ((COUNT += 1))
    RUN_DIR="$(dirname -- "$INPUT")"
    RUN_NAME="$(basename -- "$RUN_DIR")"
    RUN_START=$SECONDS

    printf "\n${BOLD}${YELLOW}[%d/%d]${RESET} ${BOLD}%s${RESET}\n" \
        "$COUNT" "$TOTAL" "$RUN_NAME"

    printf "${CYAN}  Command:${RESET} lmp -k on g 1 -sf kk -pk kokkos neigh half newton on -echo none -in in.lammps\n"

    if (
        cd -- "$RUN_DIR"
        lmp -k on g 1 -sf kk \
            -pk kokkos neigh half newton on \
            -echo none -in in.lammps \
            > lammps.stdout.log 2> lammps.stderr.log
    ); then
        ELAPSED=$((SECONDS - RUN_START))
        printf "${GREEN}${BOLD}  SUCCESS${RESET}  %s  ${GREEN}(elapsed %s)${RESET}\n" \
            "$RUN_NAME" "$(format_time "$ELAPSED")"
    else
        STATUS=$?
        ELAPSED=$((SECONDS - RUN_START))
        ((FAILED += 1))

        printf "${RED}${BOLD}  FAILED${RESET}   %s  ${RED}(exit %d, elapsed %s)${RESET}\n" \
            "$RUN_NAME" "$STATUS" "$(format_time "$ELAPSED")"
        printf "  ${YELLOW}Logs:${RESET} %s\n" "$RUN_DIR/lammps.stderr.log"
    fi
done

BATCH_ELAPSED=$((SECONDS - BATCH_START))

printf "\n${BOLD}${BLUE}========================================${RESET}\n"
printf "${BOLD} Batch summary${RESET}\n"
printf "  Total:     %d\n" "$TOTAL"
printf "  ${GREEN}Successful: %d${RESET}\n" "$((TOTAL - FAILED))"
printf "  ${RED}Failed:     %d${RESET}\n" "$FAILED"
printf "  Total time: ${BOLD}%s${RESET}\n" "$(format_time "$BATCH_ELAPSED")"
printf "${BOLD}${BLUE}========================================${RESET}\n"

if (( FAILED > 0 )); then
    exit 1
fi
