const INTERVAL = 50;

class Processor extends AudioWorkletProcessor {
  constructor() {
    super();
    this.buffer = [];
    this.sampleRate = sampleRate;
    this.targetLength = Math.floor(this.sampleRate * (INTERVAL / 1000));  }

  process(inputs) {
    const input = inputs[0][0];
    if (!input) return true;

    for (let i = 0; i < input.length; i++) {
      this.buffer.push(input[i]);
    }

    while (this.buffer.length >= this.targetLength) {
      const snippet = this.buffer.slice(0, this.targetLength);
      this.buffer = this.buffer.slice(this.targetLength);
      this.port.postMessage(snippet);
    }

    return true;
  }
}

registerProcessor('processor', Processor);
