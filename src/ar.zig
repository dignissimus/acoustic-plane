const std = @import("std");
const zpoly = @import("zpoly");

const ally = std.heap.wasm_allocator;
const RATE = 44100;
export fn burgs_method(sample: [*]f32, coefficients: [*]f32, rval: [*]f32, sample_size: usize, order: usize) i8 {
    const F = ally.alloc(f32, sample_size) catch return -1;
    const B = ally.alloc(f32, sample_size) catch return -2;
    const roots = ally.alloc(std.math.Complex(f32), order - 1) catch return -3;
    const is_good = ally.alloc(bool, order - 1) catch return -4;
    const formants = ally.alloc(f32, order - 1) catch return -5;

    defer ally.free(F);
    defer ally.free(B);
    defer ally.free(roots);
    defer ally.free(is_good);
    defer ally.free(formants);

    for (0..order) |i| {
        coefficients[i] = 0;
    }
    coefficients[0] = 1;
    for (0..sample_size) |i| {
        F[i] = sample[i];
        B[i] = sample[i];
    }

    var Dk: f32 = 0;
    for (0..sample_size) |i| {
        Dk += 2 * F[i] * F[i];
    }
    Dk -= sample[0] * sample[0] + sample[sample_size - 1] * sample[sample_size - 1];
    for (0..order - 1) |k| {
        var mu: f32 = 0;
        for (0..sample_size - k - 1) |n| {
            mu += F[n + k + 1] * B[n];
        }
        mu *= -2 / Dk;
        for (0..(k + 1) / 2 + 1) |n| {
            const t1: f32 = coefficients[n] + mu * coefficients[k + 1 - n];
            const t2: f32 = coefficients[k + 1 - n] + mu * coefficients[n];
            coefficients[n] = t1;
            coefficients[k + 1 - n] = t2;
        }
        for (0..sample_size - k - 1) |n| {
            const delta1 = mu * B[n];
            const delta2 = mu * F[n + k + 1];
            F[n + k + 1] += delta1;
            B[n] += delta2;
        }
        Dk = (1.0 - mu * mu) * Dk - F[k + 1] * F[k + 1] - B[sample_size - k - 2] * B[sample_size - k - 2];
    }

    zpoly.roots(f32, ally, coefficients[0..order], roots[0 .. order - 1]) catch return -6;
    for (roots[0 .. order - 1], is_good[0 .. order - 1]) |r, *ig| {
        ig.* = r.im > 0;
    }
    for (roots[0 .. order - 1], formants[0 .. order - 1], is_good[0 .. order - 1]) |*r, *f, ig| {
        const bandwidth = -1 / 2 * RATE / (2 * std.math.pi) * std.math.log(f32, std.math.e, std.math.complex.abs(r.*));
        const formant = std.math.complex.arg(r.*) / (2 * std.math.pi) * RATE;
        if (ig and bandwidth < 400 and formant > 90) {
            f.* = formant;
        } else {
            f.* = std.math.inf(f32);
        }
    }
    std.mem.sort(f32, formants[0 .. order - 1], {}, comptime std.sort.asc(f32));
    rval[0] = formants[0];
    rval[1] = formants[1];
    return 0;
}

export fn alloc(size: usize) ?[*]f32 {
    return if (std.heap.wasm_allocator.alloc(f32, size)) |slice| slice.ptr else |_| null;
}

export fn free(ptr: ?[*]const f32, size: usize) void {
    if (ptr) |valid_ptr| std.heap.wasm_allocator.free(valid_ptr[0..size]);
}
