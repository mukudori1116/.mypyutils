from requests_html import HTMLSession

url = "http://www.meigensyu.com/quotations/index/random"
session = HTMLSession()
res = session.get(url)
phrase = res.html.find(
    '#contents_box > div:nth-child(8) > div.text', first=True).text
origin = res.html.find(
    '.link > ul > li > a', first=True).text

print(phrase)
print(f"～{origin}～")
