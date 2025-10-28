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
} from '@fortawesome/free-solid-svg-icons'
import SvgIcon from '@mui/material/SvgIcon'

function FAIcon(props) {
  return (
    <SvgIcon
      viewBox="0 0 530 530"
      style={{
        position: 'relative',
        top: -1,
        height: '1rem',
      }}
      {...props}
    >
      <path d={props.icon.icon[4]} />
    </SvgIcon>
  )
}

function categoryIcon(category, color = '#aaa', props) {
  let iconClass = faQuestionCircle
  const iconMap = {
    P: faDna,
    D: faFlask,
    F: faAdjust,
    L: faLightbulb,
    C: faCamera,
    CF: faSliders,
    CL: faBolt,
    '%': faPercent,
  }
  if (category in iconMap) {
    iconClass = iconMap[category]
  }
  return <FAIcon icon={iconClass} htmlColor={color} {...props} />
}

export { categoryIcon, FAIcon }
