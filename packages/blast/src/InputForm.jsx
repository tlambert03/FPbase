import React from "react"
import Cookies from "js-cookie"
import Form from "react-bootstrap/Form"
import Button from "react-bootstrap/Button"

function InputForm({ onSubmit }) {
  function handleSubmit(e) {
    e.preventDefault()
    onSubmit(e.target)
  }

  return (
    <Form onSubmit={e => handleSubmit(e)}>
      <input
        type="hidden"
        name="csrfmiddlewaretoken"
        value={Cookies.get("csrftoken")}
      />
      <Form.Group controlId="queryInput">
        <Form.Label>Enter Query</Form.Label>
        <Form.Control required name="query" as="textarea" rows="4" />
        <Form.Text className="text-muted">
          Single sequence or multiple sequences in FASTA format. Accepts either
          amino acid or nucleotide sequences, but all sequences in a query be of
          the same type.
        </Form.Text>
      </Form.Group>
      <Button variant="secondary" type="submit">
        Submit
      </Button>
    </Form>
  )
}

export default InputForm
