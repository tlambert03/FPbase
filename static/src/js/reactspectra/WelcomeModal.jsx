import React, { useState, useEffect } from "react";
import Modal from "@material-ui/core/Modal";
import Typography from "@material-ui/core/Typography";
import { makeStyles } from "@material-ui/core/styles";
import Checkbox from "@material-ui/core/Checkbox";
import FormGroup from "@material-ui/core/FormGroup";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import Button from "@material-ui/core/Button";
import Box from "@material-ui/core/Box";
import List from "@material-ui/core/List";
import ListItem from "@material-ui/core/ListItem";
import ListItemIcon from "@material-ui/core/ListItemIcon";
import ListItemText from "@material-ui/core/ListItemText";
import Divider from "@material-ui/core/Divider";
import SearchIcon from "@material-ui/icons/Search";
import ChartIcon from "@material-ui/icons/InsertChart";
import FileIcon from "@material-ui/icons/GetApp";
import CachedIcon from "@material-ui/icons/Cached";

function getModalStyle() {
  const top = 40;
  const left = 42;

  return {
    top: `${top}%`,
    left: `${left}%`,
    transform: `translate(-${top}%, -${left}%)`
  };
}

const useStyles = makeStyles(theme => ({
  root: {
    width: "100%",
    maxWidth: 500
  },
  paper: {
    position: "absolute",
    width: "80%",
    backgroundColor: theme.palette.background.paper,
    boxShadow: theme.shadows[5],
    padding: 30,
    paddingTop: 25,
    outline: "none"
  },
  button: {
    height: 40,
    marginRight: 20
  },
  title: {
    color: "black",
    marginBottom: 15
  }
}));

const SearchModal = () => {
  const [modalStyle] = useState(getModalStyle);
  const classes = useStyles();
  const storageKey = "_hideFPbaseSpectraWelcome";
  const hide = localStorage.getItem(storageKey) === "true";
  const [checked, setChecked] = useState(hide);
  const [open, setOpen] = useState(!hide);

  const handleChange = e => {
    localStorage.setItem(storageKey, e.target.checked);
    setChecked(e.target.checked);
  };

  return (
    <Modal
      aria-labelledby="welcome-modal"
      aria-describedby="welcome-modal"
      open={open}
      onClose={() => setOpen(false)}
    >
      <div style={modalStyle} className={classes.paper}>
        <Typography
          variant="h6"
          id="welcome-modal-title"
          className={classes.title}
        >
          Welcome to the new FPbase Spectra Viewer!
        </Typography>

        <div className={classes.root}>
          <Typography variant="h4" gutterBottom>
            h4. Heading
          </Typography>
        </div>

        <List>
          <ListItem>
            <ListItemIcon>
              <SearchIcon />
            </ListItemIcon>
            <ListItemText primary="Quick Spectrum Entry" />
          </ListItem>
          <ListItem>
            <ListItemText inset primary="Quick Spectrum Entry" />
          </ListItem>
          <ListItem>
            <ListItemIcon>
              <ChartIcon />
            </ListItemIcon>
            <ListItemText primary="Improved graphics" />
          </ListItem>
          <ListItem>
            <ListItemIcon>
              <FileIcon />
            </ListItemIcon>
            <ListItemText primary="Exporting" />
          </ListItem>
          <ListItem>
            <ListItemIcon>
              <CachedIcon />
            </ListItemIcon>
            <ListItemText primary="State Recovery" />
          </ListItem>
        </List>

        <Box display="flex">
          <Button
            variant="contained"
            color="primary"
            className={classes.button}
            onClick={() => setOpen(false)}
          >
            Got it
          </Button>
          <FormGroup row>
            <FormControlLabel
              // prettier-ignore
              control={<Checkbox checked={checked} color="primary" onChange={handleChange} value="checked" />}
              label="Don't show this again"
            />
          </FormGroup>
        </Box>
      </div>
    </Modal>
  );
};

export default SearchModal;
