import React from "react"
import Typography from "@material-ui/core/Typography"
import NoSsr from "@material-ui/core/NoSsr"
import TextField from "@material-ui/core/TextField"
import Paper from "@material-ui/core/Paper"
import PropTypes from "prop-types"
import { makeStyles, useTheme } from "@material-ui/core/styles"
import "regenerator-runtime" // why do I need this?!?
import Select from "react-select"
import SortableWindowedSelect from "./SortableWindowedSelect"

const useStyles = makeStyles(theme => ({
  root: {
    flexGrow: 1,
  },
  input: {
    display: "flex",
    padding: 0,
    height: "auto",
    paddingBottom: 4,
    fontSize: "1.2rem",
  },
  valueContainer: {
    display: "flex",
    flexWrap: "wrap",
    flex: 1,
    alignItems: "center",
    overflow: "hidden",
  },
  noOptionsMessage: {
    padding: theme.spacing(1, 2),
  },
  singleValue: {
    fontSize: 16,
  },
  placeholder: {
    position: "absolute",
    left: 4,
    bottom: 7,
    fontSize: "1.2rem",
  },
  paper: {
    position: "absolute",
    zIndex: 1,
    marginTop: theme.spacing(1),
    left: 0,
    right: 0,
  },
  divider: {
    height: theme.spacing(2),
  },
}))

function NoOptionsMessage({ selectProps, innerProps, children }) {
  return (
    <Typography
      color="textSecondary"
      className={selectProps.classes.noOptionsMessage}
      {...innerProps}
    >
      {children}
    </Typography>
  )
}

NoOptionsMessage.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  selectProps: PropTypes.object.isRequired,
}

function inputComponent({ inputRef, ...props }) {
  return <div ref={inputRef} {...props} />
}

inputComponent.propTypes = {
  inputRef: PropTypes.oneOfType([PropTypes.func, PropTypes.object]),
}

function Control({ selectProps, innerRef, innerProps, children }) {
  return (
    <TextField
      fullWidth
      InputProps={{
        inputComponent,
        inputProps: {
          className: selectProps.classes.input,
          inputRef: innerRef,
          children,
          ...innerProps,
        },
      }}
      {...selectProps.TextFieldProps}
    />
  )
}

Control.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  innerRef: PropTypes.oneOfType([PropTypes.func, PropTypes.object]),
  selectProps: PropTypes.object.isRequired,
}

function Placeholder({ selectProps, innerProps, children }) {
  return (
    <Typography
      color="textSecondary"
      className={selectProps.classes.placeholder}
      {...innerProps}
    >
      {children}
    </Typography>
  )
}

Placeholder.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  selectProps: PropTypes.object.isRequired,
}

function SingleValue({ selectProps, innerProps, children }) {
  return (
    <Typography className={selectProps.classes.singleValue} {...innerProps}>
      {children}
    </Typography>
  )
}

function ValueContainer({ selectProps, children }) {
  return <div className={selectProps.classes.valueContainer}>{children}</div>
}

ValueContainer.propTypes = {
  children: PropTypes.node,
  selectProps: PropTypes.object.isRequired,
}

function Menu({ selectProps, innerProps, children }) {
  return (
    <Paper square className={selectProps.classes.paper} {...innerProps}>
      {children}
    </Paper>
  )
}

Menu.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  selectProps: PropTypes.object,
}

const myComponents = {
  Control,
  Menu,
  NoOptionsMessage,
  Placeholder,
  SingleValue,
  ValueContainer,
}

function MuiReactSelect({ paginate, components, ...otherprops }) {
  const classes = useStyles()
  const theme = useTheme()

  const selectStyles = {
    input: base => ({
      ...base,
      color: theme.palette.text.primary,
      "& input": {
        font: "inherit",
      },
    }),
  }

  return (
    <div className={classes.root}>
      <NoSsr>
        {paginate ? (
          <SortableWindowedSelect
            {...otherprops}
            classes={classes}
            styles={selectStyles}
            components={{ ...components, ...myComponents }}
          />
        ) : (
          <Select
            {...otherprops}
            classes={classes}
            styles={selectStyles}
            components={{ ...components, ...myComponents }}
          />
        )}
      </NoSsr>
    </div>
  )
}

MuiReactSelect.defaultProps = {
  paginate: true,
}

export default MuiReactSelect
