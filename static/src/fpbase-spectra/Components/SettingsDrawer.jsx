import React, { useEffect } from "react"
import IconButton from "@material-ui/core/IconButton"
import SettingsIcon from "@material-ui/icons/Settings"
import SwipeableDrawer from "@material-ui/core/SwipeableDrawer"
import Button from "@material-ui/core/Button"
import CloseIcon from "@material-ui/icons/Close"
import DeleteIcon from "@material-ui/icons/Cached"
import { makeStyles } from "@material-ui/core/styles"
import { useMutation } from "@apollo/react-hooks"
import gql from "graphql-tag"
import ChartOptionsForm from "./SpectraViewer/ChartOptionsForm"
import OwnerOptionsForm from "./OwnerOptionsForm"

export const useStyles = makeStyles(theme => ({
  root: {
    width: "auto",
    padding: 30,
    paddingBottom: 20,
    [theme.breakpoints.down("sm")]: {
      paddingTop: 20,
    },
  },
}))

const SettingsDrawer = () => {
  const classes = useStyles()
  const [drawerOpen, setDrawerOpen] = React.useState(false)

  const [clearForm] = useMutation(gql`
    mutation ClearForm {
      clearForm @client
    }
  `)

  useEffect(() => {
    const handleKeyDown = event => {
      // don't do anything if we're on an input
      if (
        document.activeElement &&
        document.activeElement.tagName.toUpperCase() === "INPUT"
      ) {
        return
      }
      switch (event.code) {
        case "Comma":
          event.preventDefault()
          setDrawerOpen(!drawerOpen)
          break
        default:
          break
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => {
      document.removeEventListener("keydown", handleKeyDown)
    }
  }, [drawerOpen])

  const handleClearForm = () => {
    clearForm()
    setDrawerOpen(false)
  }
  return (
    <div>
      <IconButton
        tabIndex={-1}
        onClick={() => setDrawerOpen(true)}
        edge="start"
        color="inherit"
        aria-label="Open drawer"
      >
        <SettingsIcon />
      </IconButton>
      <SwipeableDrawer
        anchor="bottom"
        open={drawerOpen}
        onClose={() => setDrawerOpen(false)}
        onOpen={() => setDrawerOpen(true)}
        BackdropProps={{ style: { backgroundColor: "rgba(0,0,0,0.2)" } }}
      >
        <div className={classes.root} role="presentation">
          <ChartOptionsForm />
          <OwnerOptionsForm />
          <Button
            color="secondary"
            onClick={() => setDrawerOpen(false)}
            style={{ marginLeft: "-25px" }}
          >
            <CloseIcon />
          </Button>
          <Button
            variant="contained"
            color="secondary"
            style={{ float: "right" }}
            onClick={handleClearForm}
          >
            <DeleteIcon />
            {' '}
&nbsp; Remove All Spectra
          </Button>
        </div>
      </SwipeableDrawer>
    </div>
  )
}

export default SettingsDrawer
