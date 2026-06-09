#!/bin/bash
# View a detector shape in the Qt GUI without the /run/reinitializeGeometry
# crash (that command is broken with ToolsSG in this Geant4 build).
#
# Sets the requested shape (and optional directional shield) as the compiled
# defaults, rebuilds, launches the GUI, and restores the previous values on exit.
#
# Usage:
#   ./view_shape.sh <shape> [shield_mm] [shield_mat]
# Examples:
#   ./view_shape.sh pP                 # 5-crystal P-pentomino, no shields
#   ./view_shape.sh tri 5 Pb           # 3-crystal petal array, 5 mm lead shields
#   ./view_shape.sh single 10 Pb       # one crystal with a 10 mm lead half-shield
set -e
SHAPE="${1:-S}"
SHIELD="${2:-0}"          # extra shield thickness in mm (0 = none)
SHIELDMAT="${3:-Al}"      # Al | Pb | W
HDR="include/DetectorConstruction.hh"
BUILD="build"

# Save current defaults.
ORIG_SHAPE=$(grep -oE 'fShape = "[^"]+"' "$HDR" | head -1 | sed 's/.*"\(.*\)"/\1/')
ORIG_THK=$(grep -oE 'fAsymCasingThk = [0-9.]+' "$HDR" | head -1 | sed 's/.*= //')
ORIG_MAT=$(grep -oE 'fAsymMaterial = "[^"]+"' "$HDR" | head -1 | sed 's/.*"\(.*\)"/\1/')

echo "Building viewer: shape=$SHAPE shield=${SHIELD}mm $SHIELDMAT"
echo "  (current defaults: $ORIG_SHAPE / ${ORIG_THK}mm / $ORIG_MAT)"

sed -i '' "s/fShape = \"$ORIG_SHAPE\"/fShape = \"$SHAPE\"/" "$HDR"
sed -i '' "s/fAsymCasingThk = $ORIG_THK/fAsymCasingThk = $SHIELD/" "$HDR"
sed -i '' "s/fAsymMaterial = \"$ORIG_MAT\"/fAsymMaterial = \"$SHIELDMAT\"/" "$HDR"

restore() {
  sed -i '' "s/fShape = \"$SHAPE\"/fShape = \"$ORIG_SHAPE\"/" "$HDR"
  sed -i '' "s/fAsymCasingThk = $SHIELD/fAsymCasingThk = $ORIG_THK/" "$HDR"
  sed -i '' "s/fAsymMaterial = \"$SHIELDMAT\"/fAsymMaterial = \"$ORIG_MAT\"/" "$HDR"
  echo "Restored defaults ($ORIG_SHAPE / ${ORIG_THK}mm / $ORIG_MAT)."
}
trap restore EXIT

cmake --build "$BUILD" >/dev/null
echo "Launching GUI (close the window to exit and restore defaults)..."
( cd "$BUILD" && ./exampleB1 )
