import React from "react"
import PropTypes from "prop-types"
import Box from "@mui/material/Box"
import ToggleButton from "@mui/material/ToggleButton"
import ToggleButtonGroup from "@mui/material/ToggleButtonGroup"
import { useQuery, useMutation } from "@apollo/client"
import { makeStyles } from "@mui/styles"
import Visibility from "@mui/icons-material/Visibility"
import { GET_ACTIVE_SPECTRA, UPDATE_ACTIVE_SPECTRA } from "../client/queries"

const useStyles = makeStyles(theme => ({
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
  const TYPE_ORDER = ["AB", "A_2P", "EX", "EM"]
  const upperA = a.subtype.toUpperCase()
  const upperB = b.subtype.toUpperCase()
  if (TYPE_ORDER.indexOf(upperA) > TYPE_ORDER.indexOf(upperB)) return 1
  return -1
}

const SubtypeSelector = React.memo(function SubtypeSelector({
  subtypes,
  skip,
}) {
  const classes = useStyles()

  const { data } = useQuery(GET_ACTIVE_SPECTRA)
  const activeSpectra = data?.activeSpectra || []

  // Create a mutable copy and add active status
  const sortedSubtypes = [...subtypes]
    .sort(subtypeSorter)
    .map(subtype => ({
      ...subtype,
      active: activeSpectra.includes(subtype.id)
    }))

  const [updateSpectra] = useMutation(UPDATE_ACTIVE_SPECTRA)
  const handleClick = e => {
    const elem = e.target.closest("button")
    const checked = !elem.classList.contains("Mui-selected")
    const variables = {}
    variables[checked ? "add" : "remove"] = [elem.value]
    updateSpectra({ variables })
  }

  // if (skip) return null

  return (
    <Box>
      <ToggleButtonGroup
        value={sortedSubtypes.filter(i => i.active).map(i => i.id)}
        size="small"
        className={classes.toggleButtonGroup}
      >
        {sortedSubtypes.map(st => (
          <ToggleButton
            key={st.id}
            value={st.id}
            onClick={handleClick}
            className={classes.toggleButton}
            tabIndex={-1}
          >
            {skip ? <Visibility /> : st.subtype.replace(/^A_/g, "")}
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
