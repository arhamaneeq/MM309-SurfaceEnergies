
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

# Resolve paths independently of the invocation directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"
RUNS_DIR="$SCRIPT_DIR/runs/bulk_equilibration"

# Permit overriding the executable without editing this script.
LMP_BIN="${LMP_BIN:-lmp}"

cd "$REPO_ROOT"

if [[ ! -d "$RUNS_DIR" ]]; then
    printf "${RED}Error:${RESET} bulk equilibration runs not found: %s\n" \
        "$RUNS_DIR" >&2
    exit 1
fi

shopt -s nullglob
INPUTS=("$RUNS_DIR"/*/in.lammps)

if (( ${#INPUTS[@]} == 0 )); then
    printf "${RED}Error:${RESET} no in.lammps files found under %s\n" \
        "$RUNS_DIR" >&2
    exit 1
fi

TOTAL=${#INPUTS[@]}
FAILED=0
SUCCEEDED=0
SKIPPED=0
COUNT=0
BATCH_START=$SECONDS

# Partition jobs BEFORE preflight. Never run -skiprun on completed jobs.
PENDING_INPUTS=()

for INPUT in "${INPUTS[@]}"; do
    RUN_DIR="$(dirname -- "$INPUT")"

    if [[ -f "$RUN_DIR/COMPLETE" ]]; then
        ((SKIPPED += 1))
    else
        PENDING_INPUTS+=("$INPUT")
    fi
done

REMAINING=${#PENDING_INPUTS[@]}

printf "\n${BOLD}${CYAN}========================================${RESET}\n"
printf "${BOLD}${CYAN} Bulk equilibration preflight${RESET}\n"
printf "  Total jobs:       %d\n" "$TOTAL"
printf "  Already complete: %d\n" "$SKIPPED"
printf "  To validate/run:  %d\n" "$REMAINING"
printf "${BOLD}${CYAN}========================================${RESET}\n"

# Validate only jobs that do not have a COMPLETE marker.
VALIDATION_FAILED=0
VALIDATION_COUNT=0

for INPUT in "${PENDING_INPUTS[@]}"; do
    RUN_DIR="$(dirname -- "$INPUT")"
    RUN_NAME="$(basename -- "$RUN_DIR")"
    VALIDATION_LOG="$RUN_DIR/validation.log"

    printf "  Checking %s ... " "$RUN_NAME"

    if (
        cd -- "$RUN_DIR"
        "$LMP_BIN" -k on g 1 -sf kk \
            -pk kokkos neigh half newton on \
            -skiprun -echo none -in in.lammps
    ) > "$VALIDATION_LOG" 2>&1; then
        printf "${GREEN}PASS${RESET}\n"
    else
        STATUS=$?
        ((VALIDATION_FAILED += 1))
        printf "${RED}FAIL${RESET} (exit %d)\n" "$STATUS"
        printf "    Log: %s\n" "$VALIDATION_LOG"
    fi

    ((VALIDATION_COUNT += 1))
done

if (( VALIDATION_FAILED > 0 )); then
    printf "\n${RED}${BOLD}Preflight failed:${RESET} %d of %d pending inputs failed validation.\n" \
        "$VALIDATION_FAILED" "$VALIDATION_COUNT"
    printf "No simulations were launched. Inspect validation.log in the affected directories.\n"
    exit 1
fi

if (( REMAINING > 0 )); then
    printf "\n${GREEN}${BOLD}All %d pending inputs passed preflight.${RESET}\n" \
        "$VALIDATION_COUNT"
else
    printf "\n${GREEN}All jobs are already marked complete; no preflight or simulations needed.${RESET}\n"
fi

# Format elapsed seconds as HH:MM:SS.
format_time() {
    local elapsed=$1
    printf '%02d:%02d:%02d' \
        "$((elapsed / 3600))" \
        "$(((elapsed % 3600) / 60))" \
        "$((elapsed % 60))"
}

printf "\n${BOLD}${BLUE}========================================${RESET}\n"
printf "${BOLD}${BLUE} Bulk equilibration batch${RESET}\n"
printf "  Total jobs:       %d\n" "$TOTAL"
printf "  Already complete: %d\n" "$SKIPPED"
printf "  To run:           %d\n" "$REMAINING"
printf "${BOLD}${BLUE}========================================${RESET}\n"

# Execute only the pending jobs that passed preflight.
for INPUT in "${PENDING_INPUTS[@]}"; do
    RUN_DIR="$(dirname -- "$INPUT")"
    RUN_NAME="$(basename -- "$RUN_DIR")"

    ((COUNT += 1))
    RUN_START=$SECONDS

    printf "\n${BOLD}${YELLOW}[%d/%d]${RESET} ${BOLD}%s${RESET}\n" \
        "$COUNT" "$REMAINING" "$RUN_NAME"

    printf "${CYAN}  Command:${RESET} %s -k on g 1 -sf kk -pk kokkos neigh half newton on -echo none -in in.lammps\n" \
        "$LMP_BIN"

    if (
        cd -- "$RUN_DIR"
        "$LMP_BIN" -k on g 1 -sf kk \
            -pk kokkos neigh half newton on \
            -echo none -in in.lammps \
            > lammps.stdout.log 2> lammps.stderr.log
    ); then
        touch "$RUN_DIR/COMPLETE"
        ((SUCCEEDED += 1))

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
printf "${BOLD} Bulk equilibration summary${RESET}\n"
printf "  Total jobs:       %d\n" "$TOTAL"
printf "  Already complete: %d\n" "$SKIPPED"
printf "  Attempted:        %d\n" "$COUNT"
printf "  ${GREEN}Succeeded:        %d${RESET}\n" "$SUCCEEDED"
printf "  ${RED}Failed:            %d${RESET}\n" "$FAILED"
printf "  Total time:       ${BOLD}%s${RESET}\n" \
    "$(format_time "$BATCH_ELAPSED")"
printf "${BOLD}${BLUE}========================================${RESET}\n"

if (( FAILED > 0 )); then
    exit 1
fi
