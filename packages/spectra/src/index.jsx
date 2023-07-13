import React, { useRef } from "react"
import "./index.css"
import { ApolloProvider } from "@apollo/react-hooks"
import App from "./App"
import { SpectraViewerContainer } from "./Components/SpectraViewer"
import initializeClient from "./client/client"
import { defaults } from "./client/resolvers"


const AppWrapper = ({ uri }) => {
  const client = useRef(initializeClient({ uri }))
  return (
    <ApolloProvider client={client.current}>
      <App />
    </ApolloProvider>
  )
}

AppWrapper.defaultProps = {
  uri: "/graphql/",
}

export default AppWrapper

const SimpleSpectraViewer = ({ uri, ids, overlaps, options, hidden }) => {
  const client = useRef(initializeClient({ uri }))

  return (
    <ApolloProvider client={client.current}>
      <Inner ids={ids} overlaps={overlaps} options={options} hidden={hidden} />
    </ApolloProvider>
  )
}

SimpleSpectraViewer.defaultProps = {
  uri: "/graphql/",
}

const Inner = ({ ids = [], overlaps = [], options, hidden = [] }) => {
  const provideOptions = Object.assign(defaults.chartOptions, options)

  return (
    <SpectraViewerContainer
      provideSpectra={ids.map(String)}
      provideOverlaps={overlaps}
      provideOptions={provideOptions}
      provideHidden={hidden.map(String)}
    />
  )
}

export { SimpleSpectraViewer }
