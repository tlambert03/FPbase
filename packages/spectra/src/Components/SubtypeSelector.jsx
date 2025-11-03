import Visibility from "@mui/icons-material/Visibility"
import Box from "@mui/material/Box"
import ToggleButton from "@mui/material/ToggleButton"
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup"
import { makeStyles } from "@mui/styles"
import PropTypes from "prop-types"
import React from "react"
import { useSpectraStore } from "../store/spectraStore"

const useStyles = makeStyles((theme) => ({
  toggleButton: {
    paddingLeft: "11px",
    paddingRight: "11px",
    height: "38px",
    // [theme.breakpoints.down(960)]: {
    //   height: "34px"
    // },
  },
  toggleButtonGroup: {
    marginLeft: "5px",
    [theme.breakpoints.down("xs")]: {
      display: "none",
    },
  },
}))

function subtypeSorter(a, b) {
  const TYPE_ORDER = ["AB", "2P", "EX", "EM"]
  const upperA = a.subtype.toUpperCase()
  const upperB = b.subtype.toUpperCase()
  if (TYPE_ORDER.indexOf(upperA) > TYPE_ORDER.indexOf(upperB)) return 1
  return -1
}

const SubtypeSelector = React.memo(function SubtypeSelector({ subtypes, skip }) {
  const classes = useStyles()

  const activeSpectra = useSpectraStore((state) => state.activeSpectra)
  const hiddenSpectra = useSpectraStore((state) => state.hiddenSpectra)
  const toggleSpectrumVisibility = useSpectraStore((state) => state.toggleSpectrumVisibility)

  // Create a mutable copy and add active status (visible = in activeSpectra AND not hidden)
  const sortedSubtypes = [...subtypes].sort(subtypeSorter).map((subtype) => ({
    ...subtype,
    active: activeSpectra.includes(subtype.id) && !hiddenSpectra.includes(subtype.id),
  }))

  const handleClick = (e) => {
    const elem = e.target.closest("button")
    const spectrumId = elem.value
    // Toggle visibility instead of adding/removing from activeSpectra
    toggleSpectrumVisibility(spectrumId)
  }

  // if (skip) return null

  return (
    <Box>
      <ToggleButtonGroup
        value={sortedSubtypes.filter((i) => i.active).map((i) => i.id)}
        size="small"
        className={classes.toggleButtonGroup}
      >
        {sortedSubtypes.map((st) => (
          <ToggleButton
            key={st.id}
            value={st.id}
            onClick={handleClick}
            className={classes.toggleButton}
            tabIndex={-1}
          >
            {skip ? <Visibility /> : st.subtype}
          </ToggleButton>
        ))}
      </ToggleButtonGroup>
    </Box>
  )
})

SubtypeSelector.propTypes = {
  subtypes: PropTypes.arrayOf(
    PropTypes.objectOf(PropTypes.oneOfType([PropTypes.string, PropTypes.bool]))
  ).isRequired,
}

export default SubtypeSelector
