#!/usr/bin/env bash
# ============================================
#  Setup Script for Linux / macOS
#  Cell Annealing Tool - Symmetric Relaxation of 2D Cellular Networks
#  Requires Python 3.14+ managed via uv
# ============================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "============================================"
echo " Cell Annealing Tool - Environment Setup"
echo "============================================"
echo ""

# ---- Step 1: Check uv, install Python 3.14.6, create .venv ----
echo "[1/4] Checking uv installation..."
if ! command -v uv &> /dev/null; then
    echo "[ERROR] uv is not found. Please install uv first:"
    echo "        https://docs.astral.sh/uv/getting-started/installation/"
    exit 1
fi
echo "[OK] uv detected:"
uv --version
echo ""
echo "      Installing Python 3.14.6 and creating .venv..."
uv python install 3.14.6
[ -d ".venv" ] || uv venv .venv --python 3.14.6
echo "[OK] .venv is ready."

# ---- Step 2: Install Python dependencies ----
echo ""
echo "[2/4] Installing Python dependencies..."
uv pip install -r requirements.txt
echo "[OK] Python dependencies installed."

# ---- Step 3: Check / Setup R ----
echo ""
echo "[3/4] Checking R environment..."

# ---- function: install R packages ----
install_r_packages() {
    local rscript_cmd="$1"
    echo ""
    echo "[4/4] Installing required R packages (conicfit, sp, shotGroups)..."
    "$rscript_cmd" -e "
        packages <- c('conicfit', 'sp', 'shotGroups')
        for (pkg in packages) {
            if (!require(pkg, character.only = TRUE)) {
                install.packages(pkg, repos = 'https://cloud.r-project.org')
            }
        }
        cat('[OK] All R packages installed.\n')
    "
    echo "[OK] All R packages installed successfully."
}

# Try using bundled R_Dist first
if [ -f "$SCRIPT_DIR/R_Dist/bin/R" ]; then
    echo "[INFO] Found bundled R_Dist."
    export R_HOME="$SCRIPT_DIR/R_Dist"
    export PATH="$SCRIPT_DIR/R_Dist/bin:$PATH"
    if [ -d "$SCRIPT_DIR/R_Dist/bin/x64" ]; then
        export PATH="$SCRIPT_DIR/R_Dist/bin/x64:$PATH"
    fi
    echo "[OK] R detected (bundled R_Dist):"
    R --version | head -n 1
    install_r_packages "Rscript"
elif command -v R &> /dev/null; then
    echo "[INFO] Using system R."
    echo "[OK] R detected:"
    R --version | head -n 1
    install_r_packages "Rscript"
else
    echo "[WARNING] R is not found (neither bundled R_Dist nor system R)."
    echo ""
    echo "  Attempting to install R automatically..."

    # Detect OS
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "[INFO] Linux detected. Installing R via package manager..."
        if command -v apt &> /dev/null; then
            echo "  Using apt (Debian/Ubuntu)..."
            sudo apt update
            sudo apt install -y r-base r-base-dev
        elif command -v dnf &> /dev/null; then
            echo "  Using dnf (Fedora)..."
            sudo dnf install -y R R-devel
        elif command -v yum &> /dev/null; then
            echo "  Using yum (CentOS/RHEL)..."
            sudo yum install -y R R-devel
        elif command -v pacman &> /dev/null; then
            echo "  Using pacman (Arch)..."
            sudo pacman -S --noconfirm r
        else
            echo "[ERROR] Unsupported Linux distribution."
            echo "  Please install R manually from: https://www.r-project.org/"
            exit 1
        fi
        echo "[OK] R installed."
        install_r_packages "Rscript"

    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo "[INFO] macOS detected."
        if command -v brew &> /dev/null; then
            echo "  Using Homebrew..."
            brew install --cask r
            echo "[OK] R installed."
            install_r_packages "Rscript"
        else
            echo "[WARNING] Homebrew not found."
            echo "  Attempting to download R from CRAN..."
            R_PKG_URL="https://cloud.r-project.org/bin/macosx/base/R-4.4.2.pkg"
            R_PKG_PATH="/tmp/R-4.4.2.pkg"
            echo "  Downloading R-4.4.2.pkg..."
            curl -L -o "$R_PKG_PATH" "$R_PKG_URL"
            echo "  Installing R (may require password)..."
            sudo installer -pkg "$R_PKG_PATH" -target /
            rm "$R_PKG_PATH"
            echo "[OK] R installed."
            install_r_packages "Rscript"
        fi
    else
        echo "[ERROR] Unsupported OS ($OSTYPE)."
        echo "  Please install R manually from: https://www.r-project.org/"
        echo "  After installing R, run the following in R console:"
        echo "      install.packages(c(\"conicfit\", \"sp\", \"shotGroups\"))"
        exit 1
    fi
fi

echo ""
echo "============================================"
echo " Setup complete! You can now run:"
echo "     uv run python energy_annealing_main.py"
echo "============================================"
