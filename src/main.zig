const std = @import("std");
const wav = @import("wav");
const zpoly = @import("zpoly");

const E: std.math.complex.Complex(f32) = .{ .re = std.math.e, .im = 0 };
const I: std.math.complex.Complex(f32) = .{ .re = 0, .im = 1 };
const PI: std.math.complex.Complex(f32) = .{ .re = std.math.pi, .im = 0 };
const RATE: f32 = 44100;

fn dft(sample: []f32, freqs: []f32, amplitudes: []f32) !void {
    for (freqs, 0..) |freq, i| {
        var total: std.math.complex.Complex(f32) = .{ .re = 0, .im = 0 };
        const sample_len: f32 = @floatFromInt(sample.len);
        for (sample, 0..) |value, j| {
            const complex_value: std.math.complex.Complex(f32) = .{ .re = value, .im = 0 };
            const factor = std.math.complex.exp(I.neg().mul(.{ .re = 2, .im = 0 }).mul(PI).mul(.{ .re = freq / RATE, .im = 0 }).mul(.{ .re = @floatFromInt(j), .im = 0 }));
            const contribution = complex_value.mul(factor).div(.{ .re = sample_len, .im = 0 });
            total = total.add(contribution);
        }
        amplitudes[i] = total.re * total.re + total.im * total.im;
    }
}

// Solve the autoregressive model using Burg's method to find the LPC coefficients
fn burgs_method(sample: []f32, coefficients: []f32, ally: std.mem.Allocator) !void {
    const F = try ally.alloc(f32, sample.len);
    const B = try ally.alloc(f32, sample.len);
    defer ally.free(F);
    defer ally.free(B);
    for (coefficients) |*coefficient| {
        coefficient.* = 0;
    }
    coefficients[0] = 1;
    for (F, B, sample) |*f, *b, amplitude| {
        f.* = amplitude;
        b.* = amplitude;
    }

    var Dk: f32 = 0;
    for (F) |f| {
        Dk += 2 * f * f;
    }
    Dk -= sample[0] * sample[0] + sample[sample.len - 1] * sample[sample.len - 1];
    for (0..coefficients.len - 1) |k| {
        var mu: f32 = 0;
        for (0..sample.len - k - 1) |n| {
            mu += F[n + k + 1] * B[n];
        }
        mu *= -2 / Dk;
        for (0..(k + 1) / 2 + 1) |n| {
            const t1: f32 = coefficients[n] + mu * coefficients[k + 1 - n];
            const t2: f32 = coefficients[k + 1 - n] + mu * coefficients[n];
            coefficients[n] = t1;
            coefficients[k + 1 - n] = t2;
        }
        for (0..sample.len - k - 1) |n| {
            const delta1 = mu * B[n];
            const delta2 = mu * F[n + k + 1];
            F[n + k + 1] += delta1;
            B[n] += delta2;
        }
        Dk = (1.0 - mu * mu) * Dk - F[k + 1] * F[k + 1] - B[sample.len - k - 2] * B[sample.len - k - 2];
    }
}

pub fn main() !void {
    const ally = std.heap.page_allocator;
    const file = try std.fs.cwd().openFile("test-cases/ee-norm-mono.wav", .{});
    defer file.close();
    var reader = std.io.bufferedReader(file.reader());
    var decoder = try wav.decoder(reader.reader());
    var data: [RATE / 20]f32 = undefined; // 1s / 20 = 50ms per frame
    var index: usize = 0;

    // lpc coefficients
    var coefficients: [50]f32 = undefined;
    var roots: [49]std.math.Complex(f32) = undefined;
    var formants: [49]f32 = undefined;
    var is_good: [49]bool = undefined;
    while (true) {
        const samples_read = try decoder.read(f32, &data);
        try burgs_method(&data, &coefficients, ally);
        index += 1;
        if (samples_read < data.len) {
            break;
        }
        try zpoly.roots(f32, ally, &coefficients, &roots);
        for (roots, &is_good) |r, *ig| {
            ig.* = r.im > 0;
        }
        for (&roots, &formants, is_good) |*r, *f, ig| {
            const bandwidth = -1 / 2 * RATE / (2 * std.math.pi) * std.math.log(f32, std.math.e, std.math.complex.abs(r.*));
            const formant = std.math.complex.arg(r.*) / (2 * std.math.pi) * RATE;
            if (ig and bandwidth < 400 and formant > 90) {
                f.* = formant;
            } else {
                f.* = std.math.inf(f32);
            }
        }
        std.mem.sort(f32, &formants, {}, comptime std.sort.asc(f32));
        std.debug.print("{}Hz, {}Hz\n", .{ formants[0], formants[1] });
    }
}

pub fn main_dft() !void {
    const ally = std.heap.page_allocator;
    const file = try std.fs.cwd().openFile("test-cases/ee-norm-mono.wav", .{});
    defer file.close();

    var reader = std.io.bufferedReader(file.reader());
    var decoder = try wav.decoder(reader.reader());

    // RATE samples for 1 second of audio
    // two channels
    var data: [RATE / 20]f32 = undefined; // 1s / 20 = 50ms per frame
    var freqs = [_]f32{ 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500 };
    const amplitudes = try ally.alloc(f32, freqs.len);
    const log_amplitudes = try ally.alloc(f32, freqs.len);
    defer ally.free(amplitudes);
    defer ally.free(log_amplitudes);

    var index: usize = 0;
    while (true) {
        index += 1;
        const samples_read = try decoder.read(f32, &data);
        try dft(data[0..samples_read], &freqs, amplitudes);

        var min_log_amplitude: f32 = 100;
        for (amplitudes, 0..) |amplitude, i| {
            const log_amplitude = std.math.log10(amplitude);
            log_amplitudes[i] = log_amplitude;
            min_log_amplitude = @min(min_log_amplitude, log_amplitude);
        }

        for (log_amplitudes) |*log_amplitude| {
            log_amplitude.* -= min_log_amplitude;
        }

        std.debug.print("t = {}ms\n", .{index * 50});
        for (freqs, log_amplitudes) |freq, amp| {
            const freq_as_int: u128 = @intFromFloat(freq);
            std.debug.print("{}: ", .{freq_as_int});
            for (0..@intFromFloat(@floor(amp))) |_| {
                std.debug.print("|", .{});
            }
            std.debug.print("\n", .{});
        }
        std.debug.print("\n Read {} samples\n", .{samples_read});

        if (samples_read < data.len) {
            break;
        }
    }
}
