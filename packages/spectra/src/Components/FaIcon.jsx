import React from "react"
import SvgIcon from "@material-ui/core/SvgIcon"
import {
  faQuestionCircle,
  faDna,
  faFlask,
  faAdjust,
  faLightbulb,
  faCamera,
  faSlidersH,
  faBolt,
  faPercent,
} from "@fortawesome/free-solid-svg-icons"

function FAIcon(props) {
  return (
    <SvgIcon
      viewBox="0 0 512 512"
      style={{
        position: "relative",
        top: -1,
        left: props.icon === faLightbulb ? -3 : -5,
        height: "1rem",
      }}
      {...props}
    >
      <path d={props.icon.icon[4]} />
    </SvgIcon>
  )
}

function categoryIcon(category, color = "#aaa", props) {
  let iconClass = faQuestionCircle
  const iconMap = {
    P: faDna,
    D: faFlask,
    F: faAdjust,
    L: faLightbulb,
    C: faCamera,
    CF: faSlidersH,
    CL: faBolt,
    "%": faPercent,
  }
  if (category in iconMap) {
    iconClass = iconMap[category]
  }
  return <FAIcon icon={iconClass} htmlColor={color} {...props} />
}

export { categoryIcon, FAIcon }
