import CloseIcon from "@mui/icons-material/Close"
import { Typography } from "@mui/material"
import IconButton from "@mui/material/IconButton"
import Paper from "@mui/material/Paper"
import Snackbar from "@mui/material/Snackbar"
import Tab from "@mui/material/Tab"
import Tabs from "@mui/material/Tabs"
import { makeStyles } from "@mui/styles"
import $ from "jquery"
import React, { useEffect, useState } from "react"
import BlastReportAlignments from "./ReportAlignments.jsx"
import BlastReportDescription from "./ReportDescription.jsx"

const useStyles = makeStyles((theme) => ({
  paperRoot: {
    width: "100%",
    marginTop: "23px",
    overflowX: "auto",
    position: "sticky",
    top: "0px",
    zIndex: "1000",
  },
  close: {
    padding: theme.spacing(0.5),
  },
}))

const NoHitsMessage = ({ open, handleClose }) => {
  const classes = useStyles()

  return (
    <div>
      <Snackbar
        anchorOrigin={{
          vertical: "bottom",
          horizontal: "left",
        }}
        open={open}
        autoHideDuration={20000}
        onClose={handleClose}
        ContentProps={{
          "aria-describedby": "message-id",
        }}
        message={
          <span id="message-id">
            <Typography key="undo" color="secondary" size="small" onClick={handleClose}>
              NO HITS
            </Typography>
            Please confirm that you are entering either amino acid sequence(s) or nucleotide
            sequence(s).
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
          </IconButton>,
        ]}
      />
    </div>
  )
}

function BlastReport({ report }) {
  const [tab, setTab] = useState(0)
  const [algnItem, setAlgnItem] = useState(null)

  function handleTabClick(_event, newValue) {
    setTab(newValue)
  }

  useEffect(() => {
    if (algnItem !== null && tab === 1) {
      $("html, body").animate(
        {
          scrollTop: $(`#dln_${algnItem}`).offset().top - 60,
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
  function closeSnackbar(_event, reason) {
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
        <Tabs value={tab} onChange={handleTabClick} indicatorColor="primary" textColor="primary">
          <Tab label="Descriptions" />
          <Tab label="Alignments" />
        </Tabs>
      </Paper>
      {tab === 0 && (
        <div>
          <BlastReportDescription report={report.report.results} onClick={handleItemClick} />
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
