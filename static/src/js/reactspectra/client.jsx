import ApolloClient from "apollo-boost";

const client = new ApolloClient({
  uri: "/graphql/"
});

window.client = client;

export default client;
