//! By convention, main.zig is where your main function lives in the case that
//! you are building an executable. If you are making a library, the convention
//! is to delete this file and start with root.zig instead.
const std = @import("std");
const wav = @import("wav");

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

fn burgs_method(sample: []f32, coefficients: []f32, ally: std.mem.Allocator) !void {
    const F = ally.alloc(f32, coefficients.size);
    const B = ally.alloc(f32, coefficients.size);
    for (&coefficients) |*coefficient| {
        coefficient.* = 0;
    }
    coefficients[0] = 1;
    for (&F, &B, sample) |*f, *b, amplitude| {
        f.* = amplitude;
        b.* = amplitude;
    }

    var Dk: f32 = -sample[0] * sample[0] - sample[sample.len - 1] * sample[sample.len - 1];
    for (F) |f| {
        Dk += 2 * f * f;
    }
    for (0..coefficients.len - 1) |k| {
        var mu: f32 = 0;
        for (0..sample.len - k - 2, B) |n, b| {
            mu += F[n + k + 1] * b;
        }
        mu *= -2 / Dk;
        for (0..(k + 1) / 2) |n| {
            const delta1: f32 = mu * coefficients[k + 1 - n];
            const delta2: f32 = mu * coefficients[n];
            coefficients[n] += delta1;
            coefficients[k + 1 - n] += delta2;
        }
        for (0..sample.len - k - 2) |n| {
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
    const file = try std.fs.cwd().openFile("test-cases/ee.wav", .{});
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

    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are belong to us.\n", .{"codebase"});

    // stdout is for the actual output of your application, for example if you
    // are implementing gzip, then only the compressed bytes should be sent to
    // stdout, not any debugging messages.
    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Run `zig build test` to run the tests.\n", .{});

    try bw.flush(); // Don't forget to flush!
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // Try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}

test "fuzz example" {
    const global = struct {
        fn testOne(input: []const u8) anyerror!void {
            // Try passing `--fuzz` to `zig build test` and see if it manages to fail this test case!
            try std.testing.expect(!std.mem.eql(u8, "canyoufindme", input));
        }
    };
    try std.testing.fuzz(global.testOne, .{});
}
