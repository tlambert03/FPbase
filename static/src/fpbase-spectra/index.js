import React, { useRef } from "react"
import "./index.css"
import App from "./App"
import { ApolloProvider } from "react-apollo-hooks"
//import { ApolloProvider } from "@apollo/react-hooks"
import initializeClient from "./client/client"

if (process.env.NODE_ENV !== "production") {
  const whyDidYouRender = require("@welldone-software/why-did-you-render")
  whyDidYouRender(React, { include: [/App.+/], logOnDifferentValues: true })
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
  uri: "/graphql/"
}

export default AppWrapper
