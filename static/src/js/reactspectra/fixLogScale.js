function fixLogScale(H) {
  // Pass error messages
  H.Axis.prototype.allowNegativeLog = true;

  // // Override conversions
  // H.Axis.prototype.log2lin = function log2lin(num) {
  //   const isNegative = num < 0;
  //   let adjustedNum = Math.abs(num);
  //   if (adjustedNum < 10) {
  //     adjustedNum += (10 - adjustedNum) / 10;
  //   }
  //   const result = Math.log(adjustedNum) / Math.LN10;
  //   return isNegative ? -result : result;
  // };

  // H.Axis.prototype.lin2log = function lin2log(num) {
  //   const isNegative = num < 0;
  //   const absNum = Math.abs(num);
  //   let result = 10 ** absNum;
  //   if (result < 10) {
  //     result = (10 * (result - 1)) / (10 - 1);
  //   }
  //   return isNegative ? -result : result;
  // };

  H.Axis.prototype.log2lin = function log2lin(num) {
    return num <= 0 ? 10 : -Math.log10(num);
  };
  
  H.Axis.prototype.lin2log = function lin2log(num) {
    return 10 ** -num;
  };
}

export default fixLogScale;
