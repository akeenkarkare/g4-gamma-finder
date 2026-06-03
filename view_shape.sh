#!/bin/bash
# View a detector shape in the Qt GUI without the /run/reinitializeGeometry
# crash (that command is broken with ToolsSG in this Geant4 build).
#
# It sets the shape as the compiled default, rebuilds, and launches the GUI.
# Restores the previous default on exit.
#
# Usage:  ./view_shape.sh S      (or square | J | T | L)
set -e
SHAPE="${1:-S}"
HDR="include/DetectorConstruction.hh"
BUILD="build"

# Save current default and swap in the requested shape.
ORIG=$(grep -oE 'fShape = "[^"]+"' "$HDR" | head -1 | sed 's/.*"\(.*\)"/\1/')
echo "Building viewer for shape=$SHAPE (current default: $ORIG)..."
sed -i '' "s/fShape = \"$ORIG\"/fShape = \"$SHAPE\"/" "$HDR"

restore() {
  sed -i '' "s/fShape = \"$SHAPE\"/fShape = \"$ORIG\"/" "$HDR"
  echo "Restored default shape to $ORIG."
}
trap restore EXIT

cmake --build "$BUILD" >/dev/null
echo "Launching GUI (close the window to exit and restore default)..."
( cd "$BUILD" && ./exampleB1 )
