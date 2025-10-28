import React, { memo } from "react"
import { makeStyles } from "@mui/styles"
import AppBar from "@mui/material/AppBar"
import Toolbar from "@mui/material/Toolbar"
import Tooltip from "@mui/material/Tooltip"
import Fab from "@mui/material/Fab"
import HelpIcon from "@mui/icons-material/Help"
import AddIcon from "@mui/icons-material/Add"

import { useMutation, useQuery } from "@apollo/client"
import FormControlLabel from "@mui/material/FormControlLabel"
import Switch from "@mui/material/Switch"
import gql from "graphql-tag"
import { IconButton } from "@mui/material"
import ShareButton from "./ShareButton"
import SettingsDrawer from "./SettingsDrawer"
import { GET_CHART_OPTIONS } from "../client/queries"
import SearchModal from "./SearchModal"

export const useStyles = makeStyles(theme => ({
  text: {
    padding: theme.spacing(2, 2, 0),
  },
  appBar: {
    top: "auto",
    bottom: 0,
    backgroundColor: "#0D4B33",
  },
  grow: {
    flexGrow: 1,
  },
  fabButton: {
    backgroundColor: "#3CA644",
    "&:hover": { backgroundColor: "#3A8B44" },
    "&:active": { backgroundColor: "#3CA644" },
    color: "white",
    position: "absolute",
    zIndex: 1,
    top: -30,
    left: 0,
    right: 0,
    margin: "0 auto",
  },
  spaceBar: {
    [theme.breakpoints.down("sm")]: {
      display: "none",
    },
    color: "rgba(255,255,255,0.2)",
    textTransform: "uppercase",
    position: "absolute",
    fontSize: "0.85rem",
    zIndex: 1,
    left: 0,
    right: 0,
    top: 30,
    margin: "0 auto",
    width: "71px",
  },
  odToggle: { paddingTop: 6, marginRight: 2 },
}))

const MyAppBar = memo(function MyAppBar({
  spectraOptions,
  clearForm,
  openHelp,
}) {
  const classes = useStyles()
  const [searchOpen, setSearchOpen] = React.useState(false)
  const handleClick = () => setSearchOpen(true)

  const { data } = useQuery(GET_CHART_OPTIONS)
  const logScale = data?.chartOptions?.logScale || false
  const [toggleLogScale] = useMutation(gql`
    mutation ToggleLogScale {
      toggleLogScale @client
    }
  `)

  return (
    <>
      <AppBar position="fixed" className={classes.appBar}>
        <Toolbar>
          <SettingsDrawer />
          <IconButton color="inherit" onClick={openHelp}>
            <HelpIcon />
          </IconButton>
          <Tooltip
            title="Click [or spacebar] for Quick Entry"
            enterDelay={700}
            leaveDelay={200}
          >
            <Fab
              tabIndex={-1}
              onClick={handleClick}
              aria-label="Add"
              className={classes.fabButton}
            >
              <AddIcon />
            </Fab>
          </Tooltip>
          <div className={classes.spaceBar}>spacebar</div>
          <div className={classes.grow} />
          <FormControlLabel
            labelPlacement="start"
            control={(
              <Switch
                tabIndex={-1}
                checked={logScale}
                onChange={toggleLogScale}
              />
)}
            label="OD"
            className={classes.odToggle}
          />
          <ShareButton />
        </Toolbar>
      </AppBar>
      {spectraOptions.length > 0 && (
        <SearchModal
          options={spectraOptions}
          open={searchOpen}
          clearForm={clearForm}
          setOpen={setSearchOpen}
        />
      )}
    </>
  )
})

export default MyAppBar
