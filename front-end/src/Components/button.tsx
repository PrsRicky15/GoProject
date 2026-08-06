import react  from "@vitejs/plugin-react-swc";
export default function  button(props:{backgroundColor:string, Value:number}) {
    const backgroundColor=props.backgroundColor;
    return <button style={{backgroundColor:props.backgroundColor}}>Value</button>
}