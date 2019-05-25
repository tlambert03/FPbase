import React, { useEffect, useContext } from "react";
import ReactDOM from "react-dom";
import SpectraSelectForm from "./SpectraSelect";
import SpectraViewer from "./SpectraViewer";
import { Store, AppContext } from "./Store";

const stateToUrl = seriesIDs => {
  return `?s=${seriesIDs.join(",")}`;
};

window.onpopstate = event => {
  console.log(event.state);
};

const App = () => {
  const { state } = useContext(AppContext);

  useEffect(() => {
    if (!state.loading) {
      const seriesIDs = state.series.map(s => s.id);
      let url = location.pathname;
      if (seriesIDs.length > 0) {
        url = `?s=${seriesIDs.join(",")}`;
      }
      console.log(url);
      window.history.pushState(seriesIDs, "", url);
    }
  }, [state]);

  return state && !state.loading ? (
    <div>
      <SpectraViewer series={state.series} />
      <SpectraSelectForm initial={state.series} initialTab={state.tab} />
    </div>
  ) : (
    <div>
      <div>loading...</div>
    </div>
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
