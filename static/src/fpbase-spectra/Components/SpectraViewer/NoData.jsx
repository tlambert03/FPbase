import React, { memo } from "react"

const NoData = memo(function NoData({ height, loading }) {
  return (
    <div
      style={{
        height: height,
        width: "100%",
        position: loading ? "relative" : "absolute",
      }}
    >
      <div
        style={{
          position: "absolute",
          top: 0.4 * height,
          zIndex: 50,
          textAlign: "center",
          fontFamily: "'Century Gothic', Roboto",
          paddingRight: "8%",
          width: "100%",
        }}
      >
        Select spectra below
        <br />
        or hit&nbsp;
        <span className="kbd">spacebar</span>
        &nbsp;for quick entry
      </div>
      <svg
        xmlns="http://www.w3.org/2000/svg"
        height={height / 1.3}
        viewBox="0 0 2947 1061"
        style={{
          fillRule: "evenodd",
          position: "absolute",
          top: 0.12 * height,
          left: "1%",
          zIndex: 49,
          margin: "auto",
          width: "100%",
        }}
      >
        <path
          style={{
            fillRule: "evenodd",
            clipRule: "evenodd",
            fill: "#000",
            opacity: 0.06,
          }}
          d="M2947,1060.88s-624.57-60.18-1103.06-405.741C1525.93,409.081,1419.45-13.191,1341.92.314c-95.24,10.978-78.62,259.289-202.86,656.871-130.38,402.445-420.5,403.7-420.5,403.7H2947Zm-1445.59,0s-201.41-.43-287.43-335.216C1110.65,323.49,1096.05,1,1047.32,1.221,977.2,1.535,902.926,557.3,647.811,769.037,308.211,1050.89,0,1060.88,0,1060.88H1501.41Z"
        />
      </svg>
    </div>
  )
})

export default NoData
