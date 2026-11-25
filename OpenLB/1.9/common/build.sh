: "${OLB_ROOT:?Set OLB_ROOT to the OpenLB 1.9r0 source tree}"
export OLB_ROOT="$(realpath "$OLB_ROOT")" OMP_NUM_THREADS=1
test -f "$OLB_ROOT/default.mk"
work="$here/runs"
mkdir -p "$work"
build=$(mktemp -d "$work/build.XXXXXX")
build_case() {
  printf 'EXAMPLE = sample\ninclude $(OLB_ROOT)/default.mk\n' > "$build/Makefile"
  nice -n 15 make -C "$build" -j1 > "$build/build.log" 2>&1 || { cat "$build/build.log" >&2; return 1; }
}
run() {
  local label=$1
  shift
  mkdir -p "$work/$label"
  (cd "$work/$label" && nice -n 15 "$build/sample" "$@" > solver.log 2>&1)
}
