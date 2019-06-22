import React from "react"
import { makeStyles } from "@material-ui/core/styles"
import AppBar from "@material-ui/core/AppBar"
import Toolbar from "@material-ui/core/Toolbar"
import Fab from "@material-ui/core/Fab"

import AddIcon from "@material-ui/icons/Add"
import SearchModal from "./SearchModal"

import { useMutation, useQuery } from "react-apollo-hooks"
import { GET_CHART_OPTIONS } from "../client/queries"
import FormControlLabel from "@material-ui/core/FormControlLabel"
import Switch from "@material-ui/core/Switch"
import gql from "graphql-tag"
import ShareButton from "./ShareButton"
import SettingsDrawer from "./SettingsDrawer"

export const useStyles = makeStyles(theme => ({
  text: {
    padding: theme.spacing(2, 2, 0)
  },
  appBar: {
    top: "auto",
    bottom: 0
  },
  grow: {
    flexGrow: 1
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
    margin: "0 auto"
  },
  spaceBar: {
    [theme.breakpoints.down("sm")]: {
      display: "none"
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
    width: "71px"
  }
}))

const MyAppBar = ({ spectraOptions, clearForm }) => {
  const classes = useStyles()
  const [searchOpen, setSearchOpen] = React.useState(false)
  const handleClick = () => setSearchOpen(true)

  const {
    data: {
      chartOptions: { logScale }
    }
  } = useQuery(GET_CHART_OPTIONS)
  const toggleLogScale = useMutation(gql`
    mutation ToggleLogScale {
      toggleLogScale @client
    }
  `)

  return (
    <>
      <AppBar position="fixed" color="primary" className={classes.appBar}>
        <Toolbar>
          <SettingsDrawer clearForm={clearForm}></SettingsDrawer>
          <Fab
            tabIndex={-1}
            onClick={handleClick}
            aria-label="Add"
            className={classes.fabButton}
          >
            <AddIcon />
          </Fab>
          <div className={classes.spaceBar}>spacebar</div>
          <div className={classes.grow} />
          <FormControlLabel
            labelPlacement="start"
            control={
              <Switch
                tabIndex={-1}
                checked={logScale}
                onChange={toggleLogScale}
                color="default"
              />
            }
            label="OD"
          />
          <ShareButton></ShareButton>
        </Toolbar>
      </AppBar>
      <SearchModal
        options={spectraOptions}
        open={searchOpen}
        clearForm={clearForm}
        setOpen={setSearchOpen}
      />
    </>
  )
}

export default MyAppBar
