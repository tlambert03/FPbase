import React, { useRef } from "react"
import "./index.css"
import { ApolloProvider } from "@apollo/react-hooks"
import App from "./App"
// import { ApolloProvider } from "@apollo/react-hooks"
import initializeClient from "./client/client"

if (process.env.NODE_ENV !== "production") {
  import("@welldone-software/why-did-you-render").then(
    ({ default: whyDidYouRender }) =>
      whyDidYouRender(React, {
        include: [],
        logOnDifferentValues: true,
        collapseGroups: true,
      })
  )
}

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
