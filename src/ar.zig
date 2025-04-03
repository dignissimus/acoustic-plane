// Solve the autoregressive model using Burg's method to find the LPC coefficients
export fn burgs_method(sample: [*]f32, coefficients: [*]f32, F: [*]f32, B: [*]f32, sample_size: u32, order: u32) void {
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
}
