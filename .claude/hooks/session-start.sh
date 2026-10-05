#!/bin/bash
# SessionStart hook: make conda available in Claude Code on the web sessions.
#
# The cloud container image is plain Ubuntu with system Python and pip -- no
# conda anywhere. This installs it so the project's own environment
# (workflow/envs/scissors.yaml) can be built, which is the only way to get the
# R side at the versions the analyses target: CRAN is blocked by the network
# policy, and Ubuntu's r-cran-* packages lag far enough to defeat the floors.
#
# Deliberately does NOT create the conda environment. That is ~400 packages and
# roughly 15 minutes, which is too long to sit in front of every session start.
# Set SCISSORS_CREATE_ENV=1 to opt in; the container state is cached after the
# hook completes, so the cost is paid once rather than per session.
#
#   bash .claude/hooks/session-start.sh                  # conda only
#   SCISSORS_CREATE_ENV=1 bash .claude/hooks/session-start.sh   # + the env

set -euo pipefail

# Local machines have their own setup; this is for the web containers only.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
    echo "session-start: not a remote session, nothing to do." >&2
    exit 0
fi

CONDA_ROOT="${CONDA_ROOT:-/opt/conda}"
# micro.mamba.pm and prefix.dev are blocked by the container network policy;
# repo.anaconda.com is reachable, so Miniconda is the route in.
INSTALLER_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

log() { echo "session-start: $*" >&2; }

# ---------------------------------------------------------------------------
# 1. Miniconda (idempotent)
# ---------------------------------------------------------------------------
if [ -x "${CONDA_ROOT}/bin/conda" ]; then
    log "conda already present at ${CONDA_ROOT} ($(${CONDA_ROOT}/bin/conda --version))"
else
    log "installing Miniconda to ${CONDA_ROOT}"
    tmp="$(mktemp -d)"
    trap 'rm -rf "${tmp}"' EXIT
    if ! curl -fsSL --max-time 600 -o "${tmp}/miniconda.sh" "${INSTALLER_URL}"; then
        log "ERROR could not download the Miniconda installer from ${INSTALLER_URL}"
        log "  The session will still work for anything that does not need conda."
        exit 0   # a missing optional tool should not block the session
    fi
    bash "${tmp}/miniconda.sh" -b -p "${CONDA_ROOT}" >/dev/null
    log "installed $(${CONDA_ROOT}/bin/conda --version)"
fi

CONDA="${CONDA_ROOT}/bin/conda"

# ---------------------------------------------------------------------------
# 2. Channels: conda-forge + bioconda only
#
# Miniconda ships Anaconda's `defaults` channels, and a current conda refuses
# to solve until their Terms of Service are accepted -- which carry
# commercial-use restrictions. workflow/envs/scissors.yaml asks for
# conda-forge and bioconda only, so drop defaults rather than accept.
# ---------------------------------------------------------------------------
"${CONDA}" config --system --remove-key channels >/dev/null 2>&1 || true
"${CONDA}" config --system --add channels bioconda
"${CONDA}" config --system --add channels conda-forge
"${CONDA}" config --system --set channel_priority strict
log "channels: $(${CONDA} config --show channels | tr -d '\n' | sed 's/channels://; s/  */ /g')"

# ---------------------------------------------------------------------------
# 3. Put conda on PATH for the session
# ---------------------------------------------------------------------------
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
    echo "export PATH=\"${CONDA_ROOT}/bin:\$PATH\"" >> "${CLAUDE_ENV_FILE}"
    log "added ${CONDA_ROOT}/bin to PATH for this session"
fi

# ---------------------------------------------------------------------------
# 4. The project environment, on request only
# ---------------------------------------------------------------------------
ENV_YAML="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}/workflow/envs/scissors.yaml"

if [ "${SCISSORS_CREATE_ENV:-0}" = "1" ]; then
    if "${CONDA}" env list | grep -qE '^scissors[[:space:]]'; then
        log "conda env 'scissors' already exists"
    elif [ -f "${ENV_YAML}" ]; then
        log "creating conda env 'scissors' (~400 packages, this takes a while)"
        "${CONDA}" env create -f "${ENV_YAML}" -n scissors >/dev/null
        log "created. activate with: conda activate scissors"
    else
        log "WARNING ${ENV_YAML} not found; skipping env creation"
    fi
else
    log "conda env not created (set SCISSORS_CREATE_ENV=1 to build it). To do it by hand:"
    log "  conda env create -f workflow/envs/scissors.yaml && conda activate scissors"
fi

log "done"
