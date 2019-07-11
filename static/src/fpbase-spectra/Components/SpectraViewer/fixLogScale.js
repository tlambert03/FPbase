function fixLogScale(H) {
  // Pass error messages
  H.Axis.prototype.allowNegativeLog = true

  H.Axis.prototype.log2lin = function log2lin(num) {
    return num <= 0 ? 10 : -Math.log10(num)
  }

  H.Axis.prototype.lin2log = function lin2log(num) {
    return 10 ** -num
  }
}

export default fixLogScale
