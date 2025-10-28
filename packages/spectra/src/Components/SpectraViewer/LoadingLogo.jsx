import React from "react"

const LoadingLogo = () => (
  <div
    style={{
      position: "relative",
      width: "100%",
      margin: "auto",
      textAlign: "center",
    }}
  >
    <span
      style={{
        top: "118px",
        left: "34%",
        position: "absolute",
        fontSize: "4.6rem",
        color: "#ccc",
      }}
    >
      loading...
    </span>
    <svg
      xmlns="http://www.w3.org/2000/svg"
      height="300"
      viewBox="0 0 2947 1061"
      style={{ fillRule: "evenodd" }}
    >
      <defs>
        <linearGradient id="lightGradient">
          <stop offset="0%" stopColor="#888">
            <animate
              attributeName="stop-color"
              values="#eee; #888"
              dur="1s"
              fill="freeze"
            />
          </stop>
          <stop offset="100%" stopColor="#eee">
            <animate
              attributeName="stop-color"
              values="#eee; #888"
              dur="1s"
              fill="freeze"
              begin="1s"
            />
          </stop>
        </linearGradient>
        <clipPath id="left-to-right">
          <rect x="0" y="0" width="0" height="1061">
            <animate
              attributeName="width"
              values="0;2947"
              dur="1.3s"
              fill="freeze"
            />
          </rect>
        </clipPath>
        <clipPath id="steady">
          <circle cx="280" cy="80" r="50" />
        </clipPath>
      </defs>
      <path
        fill="#efefef"
        d="M2947,1060.88s-624.57-60.18-1103.06-405.741C1525.93,409.081,1419.45-13.191,1341.92.314c-95.24,10.978-78.62,259.289-202.86,656.871-130.38,402.445-420.5,403.7-420.5,403.7H2947Zm-1445.59,0s-201.41-.43-287.43-335.216C1110.65,323.49,1096.05,1,1047.32,1.221,977.2,1.535,902.926,557.3,647.811,769.037,308.211,1050.89,0,1060.88,0,1060.88H1501.41Z"
      />
      <path
        fill="url(#lightGradient)"
        clipPath="url(#left-to-right)"
        d="M2947,1060.88s-624.57-60.18-1103.06-405.741C1525.93,409.081,1419.45-13.191,1341.92.314c-95.24,10.978-78.62,259.289-202.86,656.871-130.38,402.445-420.5,403.7-420.5,403.7H2947Zm-1445.59,0s-201.41-.43-287.43-335.216C1110.65,323.49,1096.05,1,1047.32,1.221,977.2,1.535,902.926,557.3,647.811,769.037,308.211,1050.89,0,1060.88,0,1060.88H1501.41Z"
      />
    </svg>
  </div>
)

export default LoadingLogo
