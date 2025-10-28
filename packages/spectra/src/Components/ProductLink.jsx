import LinkIcon from '@mui/icons-material/Link'
import IconButton from '@mui/material/IconButton'

const ProductLink = ({ current }) => {
  if (!current?.url) return null
  let ownerLink = current.url
  if (current.category === 'P') {
    ownerLink = current.url ? `/protein/${current.url}` : current.url || null
  }
  return (
    ownerLink && (
      <IconButton
        color="primary"
        aria-label="Delete"
        href={ownerLink}
        target="_blank"
        tabIndex={-1}
        style={{ padding: 6, marginLeft: 6, marginRight: 2 }}
      >
        <LinkIcon />
      </IconButton>
    )
  )
}

export default ProductLink
