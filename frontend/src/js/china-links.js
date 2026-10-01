// Links marked `data-china-only` start hidden and are shown to visitors in China.
// Pages are cached and identical for everyone, so the browser decides: by its time
// zone, or else by the country that Cloudflare's edge reports for the request.

const CHINA_TIME_ZONES = ["Asia/Shanghai", "Asia/Urumqi"]

export async function showChinaOnlyLinks() {
  const links = document.querySelectorAll("[data-china-only]")
  if (!links.length) return

  let inChina = CHINA_TIME_ZONES.includes(Intl.DateTimeFormat().resolvedOptions().timeZone)
  if (!inChina) {
    const trace = await fetch("/cdn-cgi/trace")
      .then((response) => response.text())
      .catch(() => "")
    inChina = /^loc=CN$/m.test(trace)
  }
  if (inChina) {
    for (const link of links) link.classList.remove("d-none")
  }
}
