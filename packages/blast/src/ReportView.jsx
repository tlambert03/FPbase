import React, { useState, useEffect } from "react"
import { makeStyles } from "@material-ui/core/styles"
import Paper from "@material-ui/core/Paper"
import Tabs from "@material-ui/core/Tabs"
import Tab from "@material-ui/core/Tab"
import BlastReportDescription from "./ReportDescription.jsx"
import BlastReportAlignments from "./ReportAlignments.jsx"
import Snackbar from "@material-ui/core/Snackbar"
import IconButton from "@material-ui/core/IconButton"
import CloseIcon from "@material-ui/icons/Close"
import { Typography } from "@material-ui/core"
import $ from "jquery"

const useStyles = makeStyles(theme => ({
  paperRoot: {
    width: "100%",
    marginTop: "23px",
    overflowX: "auto",
    position: "sticky",
    top: "0px",
    zIndex: "1000"
  },
  close: {
    padding: theme.spacing(0.5)
  }
}))

const NoHitsMessage = ({ open, handleClose }) => {
  const classes = useStyles()

  return (
    <div>
      <Snackbar
        anchorOrigin={{
          vertical: "bottom",
          horizontal: "left"
        }}
        open={open}
        autoHideDuration={20000}
        onClose={handleClose}
        ContentProps={{
          "aria-describedby": "message-id"
        }}
        message={
          <span id="message-id">
            <Typography
              key="undo"
              color="secondary"
              size="small"
              onClick={handleClose}
            >
              NO HITS
            </Typography>
            Please confirm that you are entering either amino acid sequence(s)
            or nucleotide sequence(s).
          </span>
        }
        action={[
          <IconButton
            key="close"
            aria-label="Close"
            color="inherit"
            className={classes.close}
            onClick={handleClose}
          >
            <CloseIcon />
          </IconButton>
        ]}
      />
    </div>
  )
}

function BlastReport({ report }) {
  const [tab, setTab] = useState(0)
  const [algnItem, setAlgnItem] = useState(null)

  function handleTabClick(event, newValue) {
    setTab(newValue)
  }

  useEffect(() => {
    if (algnItem !== null && tab === 1) {
      $("html, body").animate(
        {
          scrollTop: $("#dln_" + algnItem).offset().top - 60
        },
        300
      )
      setAlgnItem(null)
    }
  }, [algnItem, tab])

  function handleItemClick(event) {
    event.preventDefault()
    setTab(1)
    setAlgnItem(event.target.getAttribute("href"))
  }

  const classes = useStyles()

  const [open, setOpen] = React.useState(true)
  function closeSnackbar(event, reason) {
    if (reason === "clickaway") {
      return
    }
    setOpen(false)
  }

  useEffect(() => {
    console.log(report.report)
    if (report.report.results.search.hits.length < 1) {
      setOpen(true)
    }
  }, [report.report, report.report.results.search.hits])

  if (report.report.results.search.hits.length < 1)
    return <NoHitsMessage open={open} handleClose={closeSnackbar} />

  return (
    <div>
      <Paper square className={classes.paperRoot}>
        <Tabs
          value={tab}
          onChange={handleTabClick}
          indicatorColor="primary"
          textColor="primary"
        >
          <Tab label="Descriptions" />
          <Tab label="Alignments" />
        </Tabs>
      </Paper>
      {tab === 0 && (
        <div>
          <BlastReportDescription
            report={report.report.results}
            onClick={handleItemClick}
          />
        </div>
      )}
      {tab === 1 && (
        <div>
          <BlastReportAlignments report={report.report.results} />
        </div>
      )}
    </div>
  )
}

export default BlastReport
