import LinkIcon from "@mui/icons-material/Link"
import IconButton from "@mui/material/IconButton"

interface ProductLinkProps {
  current?: {
    url?: string | null
    category?: string
  } | null
}

/**
 * Renders a link button to the product page (protein or other).
 * For proteins (category "P"), links to internal protein pages.
 * For other categories, uses external URLs.
 *
 * @param props - Component props
 * @returns Link button or null if no URL available
 */
const ProductLink = ({ current }: ProductLinkProps) => {
  if (!current?.url) return null
  let ownerLink: string | null = current.url
  if (current.category === "P") {
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
