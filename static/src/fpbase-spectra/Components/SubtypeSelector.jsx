import React from "react"
import PropTypes from "prop-types"
import Box from "@material-ui/core/Box"
import ToggleButton from "@material-ui/lab/ToggleButton"
import ToggleButtonGroup from "@material-ui/lab/ToggleButtonGroup"
import { useQuery, useMutation } from "react-apollo-hooks"
import { GET_ACTIVE_SPECTRA, UPDATE_ACTIVE_SPECTRA } from "../client/queries"
import { makeStyles } from "@material-ui/core"
import Visibility from "@material-ui/icons/Visibility"

const useStyles = makeStyles(theme => ({
  toggleButton: {
    height: "38px",
    [theme.breakpoints.down("sm")]: {
      height: "34px"
    },
    [theme.breakpoints.down("xs")]: {
      display: "none"
    }
  },
  toggleButtonGroup: { marginLeft: "5px" }
}))

function subtypeSorter(a, b) {
  const TYPE_ORDER = ["AB", "A_2P", "EX", "EM"]
  const upperA = a.subtype.toUpperCase()
  const upperB = b.subtype.toUpperCase()
  if (TYPE_ORDER.indexOf(upperA) > TYPE_ORDER.indexOf(upperB)) return 1
  return -1
}

const SubtypeSelector = ({ subtypes, skip }) => {
  const classes = useStyles()

  const {
    data: { activeSpectra }
  } = useQuery(GET_ACTIVE_SPECTRA)
  subtypes.sort(subtypeSorter)
  subtypes.forEach(subtype => {
    subtype.active = activeSpectra.includes(subtype.id)
  })

  const updateSpectra = useMutation(UPDATE_ACTIVE_SPECTRA)
  const handleClick = e => {
    const elem = e.target.closest("button")
    const checked = !elem.classList.contains("Mui-selected")
    const variables = {}
    variables[checked ? "add" : "remove"] = [elem.value]
    updateSpectra({ variables })
  }

  //if (skip) return null

  return (
    <Box>
      <ToggleButtonGroup
        value={subtypes.filter(i => i.active).map(i => i.id)}
        size="small"
        className={classes.toggleButtonGroup}
      >
        {subtypes.map(st => (
          <ToggleButton
            key={st.id}
            value={st.id}
            onClick={handleClick}
            className={classes.toggleButton}
            style={{ padding: 0 }}
            tabIndex={-1}
          >
            {skip ? <Visibility></Visibility> : st.subtype.replace(/^A_/g, "")}
          </ToggleButton>
        ))}
      </ToggleButtonGroup>
    </Box>
  )
}
SubtypeSelector.propTypes = {
  subtypes: PropTypes.arrayOf(
    PropTypes.objectOf(PropTypes.oneOfType([PropTypes.string, PropTypes.bool]))
  ).isRequired
}

export default SubtypeSelector
