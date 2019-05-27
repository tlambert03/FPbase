import React, { useEffect, useContext } from "react";
import ReactDOM from "react-dom";
import Button from "@material-ui/core/Button";
import SpectraSelectForm from "./SpectraSelectForm";
import SpectraViewer from "./SpectraViewer";
import { Store, AppContext } from "./Store";
import LoadingLogo from "./LoadingLogo";

const App = () => {
  const { state, dispatch } = useContext(AppContext);

  useEffect(() => {
    // window.onpopstate = event => {
    //   dispatch({
    //     type: "UPDATE",
    //     currentSpectra: event.state.currentSpectra
    //   });
    // };
  }, []); // eslint-disable-line

  useEffect(() => {
    if (!state.loading) {
      const { currentSpectra, tab } = state;
      let url = window.location.pathname;
      if (currentSpectra.length > 0) {
        url = `?s=${currentSpectra.join(",")}&t=${tab}`;
      }
      window.history.pushState(state, "", url);
    }
  }, [state.currentSpectra, state.tab]); // eslint-disable-line

  return state && !state.loading ? (
    <div>
      <SpectraViewer />
      <SpectraSelectForm />
      <Button
        variant="contained"
        color="primary"
        className="mt-2"
        onClick={() => {
          dispatch({ type: "UPDATE", currentSpectra: [] });
        }}
      >
        reset
      </Button>
    </div>
  ) : (
    <LoadingLogo />
  );
};

const initReactSpectra = elem => {
  ReactDOM.render(
    <Store>
      <App />
    </Store>,
    document.getElementById(elem)
  );
};

export default initReactSpectra;
