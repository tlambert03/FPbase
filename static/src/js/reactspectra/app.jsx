import React, { useReducer, useEffect, useContext } from "react";
import ReactDOM from "react-dom";
import SpectraSelectForm from "./SpectraSelect";
import { SpectraViewer } from "./ChartOptions";
import { reducer } from "./reducer";
import { Store, AppContext } from "./Store";

const App = () => {
  const { state } = useContext(AppContext);

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
