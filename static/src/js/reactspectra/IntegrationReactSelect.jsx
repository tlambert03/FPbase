import React, { useState } from "react";
import Select, { components } from "react-select";
import { makeStyles, useTheme } from "@material-ui/core/styles";
import Typography from "@material-ui/core/Typography";
import NoSsr from "@material-ui/core/NoSsr";
import TextField from "@material-ui/core/TextField";
import Paper from "@material-ui/core/Paper";
import MenuItem from "@material-ui/core/MenuItem";
import PropTypes from "prop-types";
import Icon from "@material-ui/core/Icon";
import { customFilterOption } from "./util";

const useStyles = makeStyles(theme => ({
  root: {
    flexGrow: 1
  },
  input: {
    display: "flex",
    padding: 0,
    height: "auto",
    paddingBottom: 4,
    fontSize: "1.2rem"
  },
  valueContainer: {
    display: "flex",
    flexWrap: "wrap",
    flex: 1,
    alignItems: "center",
    overflow: "hidden"
  },
  noOptionsMessage: {
    padding: theme.spacing(1, 2)
  },
  singleValue: {
    fontSize: 16
  },
  placeholder: {
    position: "absolute",
    left: 4,
    bottom: 7,
    fontSize: "1.2rem"
  },
  paper: {
    position: "absolute",
    zIndex: 1,
    marginTop: theme.spacing(1),
    left: 0,
    right: 0
  },
  divider: {
    height: theme.spacing(2)
  }
}));

function NoOptionsMessage({ selectProps, innerProps, children }) {
  return (
    <Typography
      color="textSecondary"
      className={selectProps.classes.noOptionsMessage}
      {...innerProps}
    >
      {children}
    </Typography>
  );
}

NoOptionsMessage.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  selectProps: PropTypes.object.isRequired
};

function inputComponent({ inputRef, ...props }) {
  return <div ref={inputRef} {...props} />;
}

inputComponent.propTypes = {
  inputRef: PropTypes.oneOfType([PropTypes.func, PropTypes.object])
};

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
          ...innerProps
        }
      }}
      {...selectProps.TextFieldProps}
    />
  );
}

Control.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  innerRef: PropTypes.oneOfType([PropTypes.func, PropTypes.object]),
  selectProps: PropTypes.object.isRequired
};

const Option = props => {
  const myProps = { ...props };
  let iconClass = "far fa-circle";
  switch (props.data.category) {
    case "P":
      iconClass = "far fa-sun";
      break;
    case "D":
      iconClass = "fa fa-sun";
      break;
    case "F":
      iconClass = "fas fa-adjust";
      break;
    case "L":
      iconClass = "fas fa-lightbulb";
      break;
    case "C":
      iconClass = "fas fa-camera";
      break;
    default:
      iconClass = "far fa-circle";
      break;
  }

  myProps.children = (
    <>
      <Icon
        style={{ fontSize: "1rem", marginRight: 8, color: "#bbb" }}
        className={iconClass}
      />
      {myProps.children}
    </>
  );
  return <components.Option {...myProps} />;
};

// function Option({ innerRef, isFocused, isSelected, innerProps, children }) {
//   return (
//     <MenuItem
//       ref={innerRef}
//       selected={isFocused}
//       component="div"
//       style={{
//         fontWeight: isSelected ? 500 : 400
//       }}
//       {...innerProps}
//     >
//       {children}
//     </MenuItem>
//   )
// }

Option.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  innerRef: PropTypes.oneOfType([PropTypes.func, PropTypes.object]),
  isFocused: PropTypes.bool,
  isSelected: PropTypes.bool
};

function Placeholder({ selectProps, innerProps, children }) {
  return (
    <Typography
      color="textSecondary"
      className={selectProps.classes.placeholder}
      {...innerProps}
    >
      {children}
    </Typography>
  );
}

Placeholder.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  selectProps: PropTypes.object.isRequired
};

function SingleValue({ selectProps, innerProps, children }) {
  return (
    <Typography className={selectProps.classes.singleValue} {...innerProps}>
      {children}
    </Typography>
  );
}

function ValueContainer({ selectProps, children }) {
  return <div className={selectProps.classes.valueContainer}>{children}</div>;
}

ValueContainer.propTypes = {
  children: PropTypes.node,
  selectProps: PropTypes.object.isRequired
};

function Menu({ selectProps, innerProps, children }) {
  return (
    <Paper square className={selectProps.classes.paper} {...innerProps}>
      {children}
    </Paper>
  );
}

const MenuList = props => {
  const { children } = props;
  return (
    // limit to 25 results
    <components.MenuList {...props}>
      {Array.isArray(children) ? children.slice(0, 25) : children}
      {children.length > 25 ? (
        <MenuItem style={{ color: "#999" }}>results truncated...</MenuItem>
      ) : null}
    </components.MenuList>
  );
};

Menu.propTypes = {
  children: PropTypes.node,
  innerProps: PropTypes.object,
  selectProps: PropTypes.object
};

const myComponents = {
  Control,
  Menu,
  MenuList,
  NoOptionsMessage,
  Option,
  Placeholder,
  SingleValue,
  ValueContainer
};

function IntegrationReactSelect({ options, defaultValue, ...otherprops }) {
  const classes = useStyles();
  const theme = useTheme();
  const [value, setValue] = useState(defaultValue);

  const selectStyles = {
    input: base => ({
      ...base,
      color: theme.palette.text.primary,
      "& input": {
        font: "inherit"
      }
    })
  };

  return (
    <div className={classes.root}>
      <NoSsr>
        <Select
          classes={classes}
          styles={selectStyles}
          options={options}
          components={myComponents}
          filterOption={customFilterOption}
          value={value}
          onChange={setValue}
          {...otherprops}
        />
      </NoSsr>
    </div>
  );
}

export default IntegrationReactSelect;
