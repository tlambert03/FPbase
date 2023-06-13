import React from "react"
import Box from "@material-ui/core/Box"
import ToggleButton from "@material-ui/lab/ToggleButton"
import ToggleButtonGroup from "@material-ui/lab/ToggleButtonGroup"
import { makeStyles, Typography } from "@material-ui/core"
import { useMutation } from "@apollo/react-hooks"
import { UPDATE_ACTIVE_SPECTRA } from "../client/queries"

const useStyles = makeStyles(theme => ({
  toggleButton: {
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

// THIS IS NOT FINISHED
const SubtypeToggle = ({ subtypes }) => {
  const classes = useStyles()

  const [updateSpectra] = useMutation(UPDATE_ACTIVE_SPECTRA)
  const handleClick = e => {
    const elem = e.target.closest("button")
    const checked = !elem.classList.contains("Mui-selected")
    const variables = {}
    variables[checked ? "add" : "remove"] = [elem.value]
    updateSpectra({ variables })
  }

  return (
    <div style={{ width: "100%", marginTop: 6 }}>
      <Box display="flex" justifyContent="flex-end" align="center">
        <Typography
          style={{
            color: "#999",
            fontSize: "0.8rem",
            marginTop: 10,
            marginRight: 6,
          }}
        >
          TOGGLE ALL
        </Typography>
        <ToggleButtonGroup
          value={subtypes}
          size="small"
          className={classes.toggleButtonGroup}
        >
          {subtypes.map(st => (
            <ToggleButton
              key={st}
              value={st}
              onClick={handleClick}
              className={classes.toggleButton}
              style={{ padding: 0 }}
              tabIndex={-1}
            >
              {st.replace(/^A_/g, "")}
            </ToggleButton>
          ))}
        </ToggleButtonGroup>
      </Box>
    </div>
  )
}

export default SubtypeToggle
