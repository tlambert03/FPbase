import NoSsr from "@mui/material/NoSsr"
import Paper from "@mui/material/Paper"
import { useTheme } from "@mui/material/styles"
import TextField from "@mui/material/TextField"
import Typography from "@mui/material/Typography"
import { makeStyles } from "@mui/styles"
import Select from "react-select"
import SortableWindowedSelect from "./SortableWindowedSelect"

const useStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
  },
  input: {
    display: "flex",
    padding: "0 4px 4px 4px",
    height: "auto",
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
    left: 8,
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

function inputComponent({ inputRef, ...props }) {
  return <div ref={inputRef} {...props} />
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

function Placeholder({ selectProps, innerProps, children }) {
  return (
    <Typography color="textSecondary" className={selectProps.classes.placeholder} {...innerProps}>
      {children}
    </Typography>
  )
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

function Menu({ selectProps, innerProps, children }) {
  return (
    <Paper square className={selectProps.classes.paper} {...innerProps}>
      {children}
    </Paper>
  )
}

const myComponents = {
  Control,
  Menu,
  NoOptionsMessage,
  Placeholder,
  SingleValue,
  ValueContainer,
}

function MuiReactSelect({ paginate = true, components, ...otherprops }) {
  const classes = useStyles()
  const theme = useTheme()

  const selectStyles = {
    input: (base) => ({
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

export default MuiReactSelect
