
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

HARD_RERUN=false

usage() {
    cat <<EOF
Usage: $0 [--hard-rerun] [--help]

Options:
  --hard-rerun  Clear generated pipeline outputs and rerun everything.
  --help        Show this help message.

Default behaviour preserves existing outputs and COMPLETE markers.
EOF
}

while (($# > 0)); do
    case "$1" in
        --hard-rerun)
            HARD_RERUN=true
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            printf "${RED}Error:${RESET} unknown argument: %s\n" "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done


# Resolve paths independently of the invocation directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd)"

cd "$REPO_ROOT"

clear_outputs() {
    printf "\n${RED}${BOLD}WARNING: HARD RERUN${RESET}\n"
    printf "This will delete generated structures, results, and simulation runs.\n"
    printf "Repository: %s\n\n" "$REPO_ROOT"

    read -r -p "Type DELETE to confirm: " CONFIRM
    if [[ "$CONFIRM" != "DELETE" ]]; then
        printf "Hard rerun cancelled. No outputs were deleted.\n"
        exit 1
    fi

    local paths=(
        "$SCRIPT_DIR/structures/bulk"
        "$SCRIPT_DIR/structures/bulk_temperature"
        "$SCRIPT_DIR/structures/slabs"
        "$SCRIPT_DIR/results"
        "$SCRIPT_DIR/runs/bulk_equilibration"
        "$SCRIPT_DIR/runs/slab_calculations"
    )

    printf "\n${YELLOW}Clearing generated outputs...${RESET}\n"

    for path in "${paths[@]}"; do
        if [[ -e "$path" ]]; then
            printf "  Removing %s\n" "$path"
            rm -rf -- "$path"
        fi
    done

    printf "${GREEN}Output cleanup complete.${RESET}\n"
}

# Preserve any executable override supplied by the caller.
export LMP_BIN="${LMP_BIN:-lmp}"

if [[ "$HARD_RERUN" == true ]]; then
    clear_outputs
else
    printf "\n${CYAN}Preserving existing outputs and completion markers.${RESET}\n"
fi

PIPELINE_START=$SECONDS
CURRENT_STAGE="initialization"

format_time() {
    local elapsed=$1
    printf '%02d:%02d:%02d' \
        "$((elapsed / 3600))" \
        "$(((elapsed % 3600) / 60))" \
        "$((elapsed % 60))"
}

run_stage() {
    local stage_name="$1"
    shift

    local stage_start=$SECONDS
    CURRENT_STAGE="$stage_name"

    printf "\n${BOLD}${BLUE}========================================${RESET}\n"
    printf "${BOLD}${BLUE} Stage: %s${RESET}\n" "$stage_name"
    printf "${BOLD}${BLUE}========================================${RESET}\n"
    printf "${CYAN}Command:${RESET}"
    printf ' %q' "$@"
    printf '\n\n'

    if "$@"; then
        local elapsed=$((SECONDS - stage_start))
        printf "\n${GREEN}${BOLD}Stage completed:${RESET} %s " "$stage_name"
        printf "${GREEN}(elapsed %s)${RESET}\n" "$(format_time "$elapsed")"
    else
        local status=$?
        local elapsed=$((SECONDS - stage_start))

        printf "\n${RED}${BOLD}PIPELINE FAILED${RESET}\n"
        printf "  Stage:  %s\n" "$stage_name"
        printf "  Exit:   %d\n" "$status"
        printf "  Elapsed: %s\n" "$(format_time "$elapsed")"
        printf "  Later stages will not run.\n"
        exit "$status"
    fi
}

require_complete_bulk_runs() {
    local runs_dir="$SCRIPT_DIR/runs/bulk_equilibration"
    local input
    local total=0
    local complete=0
    local incomplete=0

    shopt -s nullglob
    local inputs=("$runs_dir"/*/in.lammps)

    if (( ${#inputs[@]} == 0 )); then
        printf "${RED}Error:${RESET} no bulk equilibration inputs found in %s\n" \
            "$runs_dir" >&2
        return 1
    fi

    for input in "${inputs[@]}"; do
        total=$((total + 1))

        if [[ -f "$(dirname -- "$input")/COMPLETE" ]]; then
            complete=$((complete + 1))
        else
            incomplete=$((incomplete + 1))
            printf "  ${YELLOW}Incomplete:${RESET} %s\n" \
                "$(basename -- "$(dirname -- "$input")")"
        fi
    done

    printf "\nBulk completion check: %d/%d jobs complete.\n" \
        "$complete" "$total"

    if (( incomplete > 0 )); then
        printf "${RED}Error:${RESET} bulk extraction and slab generation cannot proceed.\n" \
            >&2
        printf "Finish or repair the incomplete bulk jobs first.\n" >&2
        return 1
    fi
}

printf "${BOLD}${CYAN}"
printf "========================================\n"
printf " MM309 Surface Energy Pipeline\n"
printf "========================================\n"
printf "${RESET}"
printf "Repository: %s\n" "$REPO_ROOT"
printf "LAMMPS:     %s\n" "$LMP_BIN"

# 1. Generate bulk configurations.
run_stage "Generate bulk structures" \
    python -m studies.main.bulk_structs

# 2. Generate bulk equilibration LAMMPS inputs.
run_stage "Generate bulk equilibration inputs" \
    python -m studies.main.bulk_lammps

# 3. Validate and run all bulk simulations.
run_stage "Run bulk equilibration" \
    bash "$SCRIPT_DIR/bulk_run.sh"

# 4. Ensure every bulk simulation finished before extraction.
run_stage "Verify bulk completion" \
    require_complete_bulk_runs

# 5. Extract bulk lattice parameters and reference energies.
run_stage "Extract bulk results" \
    python -m studies.main.bulk_extract

# 6. Build temperature-specific parent structures and surface slabs.
run_stage "Generate slab structures" \
    python -m studies.main.slab_struct

# 7. Generate slab LAMMPS input decks.
run_stage "Generate slab calculation inputs" \
    python -m studies.main.slab_lammps

# 8. Validate and execute slab calculations.
run_stage "Run slab calculations" \
    bash "$SCRIPT_DIR/slab_run.sh"

# Final summary.
PIPELINE_ELAPSED=$((SECONDS - PIPELINE_START))

printf "\n${BOLD}${GREEN}========================================${RESET}\n"
printf "${BOLD}${GREEN} PIPELINE COMPLETED SUCCESSFULLY${RESET}\n"
printf "${BOLD}${GREEN}========================================${RESET}\n"
printf "Total elapsed time: %s\n" "$(format_time "$PIPELINE_ELAPSED")"
printf "Bulk results:      %s/results/\n" "$SCRIPT_DIR"
printf "Slab results:      %s/runs/slab_calculations/\n" "$SCRIPT_DIR"
