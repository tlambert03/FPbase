import SvgIcon from "@mui/material/SvgIcon"
import {
  faAdjust,
  faBolt,
  faCamera,
  faDna,
  faFlask,
  faLightbulb,
  faPercent,
  faQuestionCircle,
  faSliders,
} from "../icons"

function FAIcon(props) {
  const { icon, ...rest } = props
  const viewBox = `0 0 ${icon.width} ${icon.height}`

  return (
    <SvgIcon
      viewBox={viewBox}
      style={{
        position: "relative",
        top: -1,
        height: "1rem",
      }}
      {...rest}
    >
      <path d={icon.path} />
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
    CF: faSliders,
    CL: faBolt,
    "%": faPercent,
  }
  if (category in iconMap) {
    iconClass = iconMap[category]
  }
  return <FAIcon icon={iconClass} htmlColor={color} {...props} />
}

export { categoryIcon, FAIcon }
