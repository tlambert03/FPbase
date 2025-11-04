import SvgIcon, { type SvgIconProps } from "@mui/material/SvgIcon"
import type { IconData } from "../icons"
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

interface FAIconProps extends SvgIconProps {
  icon: IconData
}

/**
 * Custom FontAwesome-style icon component using MUI's SvgIcon.
 *
 * @param props - Component props including icon data and MUI SvgIcon props
 */
function FAIcon({ icon, ...rest }: FAIconProps) {
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

/**
 * Spectrum category identifiers used throughout the application.
 */
export type SpectrumCategory = "P" | "D" | "F" | "L" | "C" | "CF" | "CL" | "%"

/**
 * Returns an icon component corresponding to a spectrum category.
 *
 * @param category - The spectrum category code
 * @param color - Icon color (default: #aaa)
 * @param props - Additional SvgIcon props
 * @returns FAIcon component with appropriate icon for the category
 */
function categoryIcon(category: string, color = "#aaa", props?: SvgIconProps) {
  let iconClass: IconData = faQuestionCircle
  const iconMap: Record<SpectrumCategory, IconData> = {
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
    iconClass = iconMap[category as SpectrumCategory]
  }
  return <FAIcon icon={iconClass} htmlColor={color} {...props} />
}

export { categoryIcon, FAIcon }
