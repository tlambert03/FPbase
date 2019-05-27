import React, { useEffect, useContext } from "react";
import ReactDOM from "react-dom";
import { ThemeProvider } from "@material-ui/styles";
import SpectraSelectForm from "./SpectraSelectForm";
import SpectraViewer from "./SpectraViewer";
import { Store, AppContext } from "./Store";
import LoadingLogo from "./LoadingLogo";
import SearchModal from "./SearchModal";
import WelcomeModal from "./WelcomeModal";
import theme from "./theme";

const App = () => {
  const { state, dispatch } = useContext(AppContext);

  useEffect(() => {
    // window.onpopstate = event => {
    //   dispatch({
    //     type: "UPDATE",
    //     currentSpectra: event.state.currentSpectra
    //   });
    // };
  }, []) // eslint-disable-line

  useEffect(() => {
    if (!state.loading && state.pushState) {
      const { currentSpectra, tab } = state;
      let url = window.location.pathname;
      if (currentSpectra.length > 0) {
        url = `?s=${currentSpectra.join(",")}&t=${tab}`;
      }
      window.history.pushState(state, "", url);
    }
  }, [state.currentSpectra, state.tab]) // eslint-disable-line

  return state && !state.loading ? (
    <div>
      <SpectraViewer />
      <SpectraSelectForm />
      <SearchModal options={Object.values(state.owners)} />
      <WelcomeModal />
    </div>
  ) : (
    <LoadingLogo />
  );
};

const initReactSpectra = elem => {
  ReactDOM.render(
    <Store>
      <ThemeProvider theme={theme}>
        <App />
      </ThemeProvider>
    </Store>,
    document.getElementById(elem)
  );
};

export default initReactSpectra;
