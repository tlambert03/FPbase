import { createMuiTheme } from "@material-ui/core/styles";

const theme = createMuiTheme({
  //   palette: {
  //     primary: purple,
  //     secondary: green
  //   },
  //   status: {
  //     danger: "orange"
  //   },
  breakpoints: {
    values: {
      xs: 0,
      sm: 768,
      md: 992,
      lg: 1280,
      sl: 1920
    }
  }
});

export default theme;
