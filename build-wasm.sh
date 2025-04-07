# zig build-exe -fno-strip src/ar.zig -target wasm32-freestanding -fno-entry --export=burgs_method
zig build -Dtarget=wasm32-freestanding
